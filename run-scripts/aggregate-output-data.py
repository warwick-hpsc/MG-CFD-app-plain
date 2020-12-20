import os
import pandas as pd
import numpy as np
import fnmatch
import argparse
import re
from math import floor, log10

import sys
pyv = sys.version_info[0]
if pyv == 2:
    from sets import Set
elif pyv == 3:
    # Create alias:
    Set = set

script_dirpath = os.path.join(os.getcwd(), os.path.dirname(__file__))
mg_cfd_dirpath = os.path.join(script_dirpath, "../")

parser = argparse.ArgumentParser()
parser.add_argument('--output-dirpath', required=True, help="Dirpath to generated processed data")
parser.add_argument('--data-dirpaths', nargs='+', default=[], required=True, help="Dirpath(s) to output data of MG-CFD runs")
parser.add_argument('--assembly-dirpath', required=False, help="Dirpath to AssemblyLoopExtractor")
args = parser.parse_args()
assembly_analyser_dirpath = args.assembly_dirpath
mg_cfd_output_dirpaths = args.data_dirpaths

prepared_output_dirpath = args.output_dirpath

if not assembly_analyser_dirpath is None:
    import imp
    imp.load_source('assembly_analysis', os.path.join(assembly_analyser_dirpath, "assembly_analysis.py"))
    from assembly_analysis import *
    imp.load_source('utils', os.path.join(assembly_analyser_dirpath, "utils.py"))
    from utils import *

compile_info = {}

essential_colnames = []
# essential_colnames += ["CPU", "CC", "CC version", "Instruction set"]
# essential_colnames += ["Event"]
# essential_colnames += ["SIMD failed", "SIMD conflict avoidance strategy"]

key_loops = ["compute_flux_edge", "unstructured_stream", "compute_stream"]

possible_data_colnames = ["Time", "Count", "NumIters", "Value"]

def grep(text, filepath):
    found = False
    f = open(filepath, "r")
    for line in f:
        if text in line:
            found = True
            break
    return found
    

def clean_pd_read_csv(filepath):
    df = pd.read_csv(filepath, keep_default_na=False)
    for v in df.columns.values:
        if v.startswith("Unnamed"):
            df = df.drop(v, axis=1)
        # elif v in ["Total", "ThreadNum", "CpuId"]:
        # elif v in ["Total", "ThreadNum"]:
        elif v in ["Total", "CpuId"]:
            df = df.drop(v, axis=1)
    return df

def safe_pd_append(top, bottom):
    top_names = top.columns.values
    bot_names = bottom.columns.values
    if Set(top_names) != Set(bot_names):
        missing_from_top = Set(bot_names).difference(top_names)
        if len(missing_from_top) > 0:
            raise Exception("Top DF missing these columns: " + missing_from_top.__str__())
        missing_from_bot = Set(top_names).difference(bot_names)
        if len(missing_from_bot) > 0:
            raise Exception("Bottom DF missing these columns: " + missing_from_bot.__str__())
    ## Ensure 'bottom' columns have same ordering:
    bottom = bottom[top_names]
    return top.append(bottom).reset_index(drop=True)

def safe_pd_filter(df, field, value):
    if not field in df.columns.values:
        print("WARNING: field '{0}' not in df".format(field))
        return df
    if isinstance(value, list):
        if len(value) == 0:
            raise Exception("No values to filter on")
        df = df[df[field].isin(value)].copy()
    else:
        df = df[df[field]==value]
        df = df.drop(field, axis=1)
    nrows = df.shape[0]
    if nrows == 0:
        raise Exception("No rows left after filter: '{0}' == '{1}'".format(field, value))
    return df

def safe_pd_merge(lhs, rhs, v=None):
    if v is None:
        d = lhs.merge(rhs)
        if d.shape[0] == 0 and (lhs.shape[0] > 0 and rhs.shape[0]):
            raise Exception("Merge produced empty DataFrame")
    else:
        shape_before = lhs.shape
        d = lhs.merge(rhs, validate=v)
        if v == "many_to_one":
            if d.shape[0] < shape_before[0]:
                raise Exception("Rows lost from LHS during merge: {0} before, {1} after".format(shape_before[0], d.shape[0]))
            elif d.shape[0] > shape_before[0]:
                raise Exception("Rows added to LHS during merge: {0} before, {1} after".format(shape_before[0], d.shape[0]))
        else:
            if d.shape[0] == 0 and (lhs.shape[0] > 0 and rhs.shape[0]):
                raise Exception("Merge produced empty DataFrame")

    return d

def safe_frame_divide(top, bottom):
    data_colnames = get_data_colnames(bottom)
    renames = {}
    for x in data_colnames:
        renames[x] = x+"_div"
    b = bottom.rename(index=str, columns=renames)
    d = top.merge(b, validate="one_to_one")
    for cn in data_colnames:
        f = d[cn+"_div"] != 0.0
        d.loc[f,cn] /= d.loc[f,cn+"_div"]
        f = d[cn+"_div"] == 0.0
        d.loc[f,cn] = 0.0
        d = d.drop(cn+"_div", axis=1)
    return d

def get_data_colnames(df):
    data_colnames = []
    for v in df.columns.values:
        if v.startswith("insn.") or v.startswith("eu."):
            data_colnames.append(v)
    data_colnames += [c for c in possible_data_colnames if c in df.columns.values]
    return data_colnames

def get_job_id_colnames(df):
    data_colnames = get_data_colnames(df)
    job_id_colnames = [c for c in df.columns.values if not c in data_colnames]

    if len(job_id_colnames) == 0:
        print("ERROR: get_job_id_colnames() failed to find any.")
        sys.exit(-1)
    return job_id_colnames

def load_attributes(output_dirpath):
    att_df = clean_pd_read_csv(os.path.join(output_dirpath, "Attributes.csv"))
    att_dict = {}
    for i in range(att_df.shape[0]):
        row = att_df.iloc[i]
        k = row["Attribute"]
        v = row["Value"]
        try:
            v = int(v)
        except:
            try:
                v = float(v)
            except:
                pass
        att_dict[k] = v
    return att_dict

def calc_ins_per_iter(output_dirpath, kernel):
    if "compute_flux_edge" in kernel:
        kernel = "compute_flux_edge"
    elif "unstructured_stream" in kernel:
        kernel = "unstructured_stream"
    elif "unstructured_compute" in kernel:
        kernel = "unstructured_compute"
    elif "compute_stream" in kernel:
        kernel = "compute_stream"
    ins_per_iter = -1

    papi_filepath = os.path.join(output_dirpath, "PAPI.csv")
    loop_num_iters_filepath = os.path.join(output_dirpath, "LoopNumIters.csv")

    if os.path.isfile(papi_filepath) and os.path.isfile(loop_num_iters_filepath):
        iters = clean_pd_read_csv(loop_num_iters_filepath)
        iters = iters[iters["Loop"]==kernel].drop("Loop", axis=1)
        simd_failed = did_simd_fail(output_dirpath)
        if simd_failed:
            simd_len = 1
        else:
            simd_len = infer_field(output_dirpath, "SIMD len")

        # 'num_iters' is counting non-vectorised iterations. If code 
        # was vectorised then fewer loop iterations will have been 
        # performed. Make that adjustment:
        iters["NumIters"] = iters["NumIters"] / simd_len

        papi = clean_pd_read_csv(papi_filepath)
        papi = papi[papi["Loop"]==kernel].drop("Loop", axis=1)
        if "PAPI_TOT_INS" in papi["Event"].unique():
            papi = papi[papi["Event"]=="PAPI_TOT_INS"]
            papi = papi.drop("Event", axis=1)

            job_id_colnames = get_job_id_colnames(papi)
            data_colnames = [c for c in papi.columns.values if not c in job_id_colnames]

            papi_and_iters = safe_pd_merge(papi, iters)
            papi_and_iters = papi_and_iters[papi_and_iters["MG level"]==0]
            papi_and_iters["InsPerIter"] = papi_and_iters["Count"] / papi_and_iters["NumIters"]
            ins_per_iter = papi_and_iters.loc[0,"InsPerIter"]

    if ins_per_iter == -1:
        print("WARNING: Failed to calculate ins-per-iter for loop '{0}'".format(kernel))

    return ins_per_iter

def infer_field(output_dirpath, field):
    value = None

    att_filepath = os.path.join(output_dirpath, "Attributes.csv")
    if not os.path.isfile(att_filepath):
        raise Exception("File not found: " + att_filepath)
    attributes = load_attributes(output_dirpath)
    if field in attributes:
        return attributes[field]

    if value is None:
        raise Exception("Cannot infer field '{0}'".format(field))
    
    return value

def did_simd_fail(output_dirpath, compiler=None):
    if infer_field(output_dirpath, "SIMD") == "N":
        return False

    failed = False
    log_filepath =  os.path.join(output_dirpath, "objects", "flux_vecloops.o.log")

    if compiler is None:
        compiler = infer_field(output_dirpath, "CC")

    if compiler == "clang":
        f = open(log_filepath, "r")
        simd_fail_found = False
        simd_success_found = False
        flux_loop_simd_fail = False
        # for line in f:
        for i,line in enumerate(f):
            if "loop not vectorized" in line:
                simd_fail_found = True
                simd_success_found = False
                ## Source code snippet should be on next line, confirming which loop
                continue
            if "vectorized loop" in line:
                simd_success_found = True
                simd_fail_found = False
                ## Source code snippet should be on next line, confirming which loop
                continue

            if "for (long i=flux_loop_start" in line:
                if simd_fail_found:
                    flux_loop_simd_fail = True
                    break
                if simd_success_found:
                    break
            simd_fail_found = False
            simd_success_found = False

        if flux_loop_simd_fail:
            failed = True

    elif compiler == "gnu":
        f = open(log_filepath, "r")
        for line in f:
            if re.match("src/Kernels/flux_kernel.elemfunc.c.*not vectorized", line):
                failed = True
                break

    return failed

def analyse_object_files():
    print("Analysing object files")

    dirpaths = mg_cfd_output_dirpaths

    for dp in dirpaths:
        for o in os.listdir(dp):
            output_dirpath = os.path.join(dp, o)
            if not os.path.isdir(output_dirpath):
                continue

            ic_filepath = os.path.join(output_dirpath, "instruction-counts.csv")
            ic_cat_filepath = os.path.join(output_dirpath, "instruction-counts-categorised.csv")
            if os.path.isfile(ic_filepath) and os.path.isfile(ic_cat_filepath):
                continue

            # print("Processing: " + output_dirpath)

            loop_to_object = {}

            run_filename_candidates = [x+".batch" for x in ["slurm", "pbs", "moab", "lsf"]]
            run_filename_candidates.append("run.sh")
            for rhc in run_filename_candidates:
                if os.path.isfile(os.path.join(output_dirpath, rhc)):
                    run_filename = rhc
                    break
            if run_filename == "":
                raise IOError("Cannot find a run script for run: " + output_dirpath)

            compile_info["compiler"] = infer_field(output_dirpath, "CC")
            compile_info["SIMD len"] = infer_field(output_dirpath, "SIMD len")
            simd = infer_field(output_dirpath, "SIMD")
            if simd == "N":
                compile_info["SIMD len"] = 1
            simd_failed = did_simd_fail(output_dirpath, compile_info["compiler"])
            compile_info["SIMD failed"] = simd_failed
            avx512cd_required = compile_info["compiler"] == "intel" and \
                                simd == "Y" and \
                                infer_field(output_dirpath, "Instruction set") == "AVX512" and \
                                infer_field(output_dirpath, "SIMD conflict avoidance strategy") == "None"
            compile_info["SIMD CA scheme"] = infer_field(output_dirpath, "SIMD conflict avoidance strategy")
            compile_info["scatter loop present"] = "manual" in compile_info["SIMD CA scheme"].lower() and "scatter" in compile_info["SIMD CA scheme"].lower()
            compile_info["gather loop present"]  = "manual" in compile_info["SIMD CA scheme"].lower() and "gather"  in compile_info["SIMD CA scheme"].lower()

            if simd == "Y":
                loop_to_object["compute_flux_edge_vecloop"] = "flux_vecloops.o"
                loop_to_object["unstructured_stream_vecloop"] = "unstructured_stream_vecloop.o"
                loop_to_object["unstructured_compute_vecloop"] = "unstructured_compute_vecloop.o"
                loop_to_object["compute_stream_vecloop"] = "compute_stream_vecloop.o"
            else:
                loop_to_object["compute_flux_edge_loop"] = "flux_loops.o"
                loop_to_object["unstructured_stream_loop"] = "unstructured_stream_loop.o"
                loop_to_object["unstructured_compute_loop"] = "unstructured_compute_loop.o"
                loop_to_object["compute_stream_loop"] = "compute_stream_loop.o"

            if "manual" in compile_info["SIMD CA scheme"].lower():
                if compile_info["compiler"] == "clang":
                    compile_info["manual CA block width"] = (compile_info["SIMD len"]*3)+2
                else:
                    compile_info["manual CA block width"] = compile_info["SIMD len"]

            loops_tally_df = None
            for k in loop_to_object.keys():
                if "compute_flux_edge" in k:
                    k_pretty = "compute_flux_edge"
                elif "unstructured_stream" in k:
                    k_pretty = "unstructured_stream"
                elif "unstructured_compute" in k:
                    k_pretty = "unstructured_compute"
                elif "compute_stream" in k:
                    k_pretty = "compute_stream"

                ins_per_iter = calc_ins_per_iter(output_dirpath, k)

                obj_filepath = os.path.join(output_dirpath, "objects", loop_to_object[k])
                try:
                    asm_loop_filepath = extract_loop_kernel_from_obj(obj_filepath, compile_info, ins_per_iter, k, avx512cd_required, 10)
                except:
                    raise
                    # continue
                loop_tally = count_loop_instructions(asm_loop_filepath)

                df_data = [ [k_pretty, k, float(v)] for k,v in loop_tally.items()]
                loop_tally_df = pd.DataFrame(df_data, columns=["Loop", "Instruction", "Count"])
                if loops_tally_df is None:
                    loops_tally_df = loop_tally_df
                else:
                    loops_tally_df = loops_tally_df.append(loop_tally_df)

            if "LOAD_SPILLS" in loops_tally_df["Instruction"].values:
                if "unstructured_stream" in loops_tally_df["Loop"].values:
                    ## This loop does not actually have any register spilling, so 
                    ## any that were identified by 'assembly-loop-extractor' are false positives. 
                    ## Assume that the same number of false positives are identified in the 'flux' kernel, 
                    ## so remove those:
                    f = np.logical_and(loops_tally_df["Loop"]=="unstructured_stream", loops_tally_df["Instruction"]=="LOAD_SPILLS")
                    claimed_unstructured_stream_spills = loops_tally_df.loc[f,"Count"].iloc[0]
                    f = loops_tally_df["Instruction"]=="LOAD_SPILLS"
                    loops_tally_df.loc[f,"Count"] -= claimed_unstructured_stream_spills

            loops_tally_df.to_csv(ic_filepath, index=False)

            target_is_aarch64 = False
            if simd == "Y":
                raw_asm_filepath = os.path.join(output_dirpath, "objects", "flux_vecloops.o.raw-asm")
            else:
                raw_asm_filepath = os.path.join(output_dirpath, "objects", "flux_loops.o.raw-asm")
            for line in open(raw_asm_filepath):
              if "aarch64" in line:
                target_is_aarch64 = True
                break

            loops_tally_categorised_df = None
            loops = Set(loops_tally_df["Loop"])
            for l in loops:
                df = loops_tally_df[loops_tally_df["Loop"]==l]
                tally = {}
                for i in range(df.shape[0]):
                    row = df.iloc[i]
                    k = row["Instruction"]
                    v = row["Count"]
                    tally[k] = v

                if target_is_aarch64:
                    tally_categorised = categorise_aggregated_instructions_tally_dict(tally, is_aarch64=True)
                else:
                    tally_categorised = categorise_aggregated_instructions_tally_dict(tally, is_intel64=True)

                df_data = [ [l, k, float(v)] for k,v in tally_categorised.items()]
                df = pd.DataFrame(df_data, columns=["Loop", "Instruction", "Count"])
                if loops_tally_categorised_df is None:
                    loops_tally_categorised_df = df
                else:
                    loops_tally_categorised_df = loops_tally_categorised_df.append(df)

            loops_tally_categorised_df.to_csv(ic_cat_filepath, index=False)

def collate_csvs():
    cats = ["Attributes", "Times", "PAPI", "instruction-counts", "instruction-counts-categorised", "LoopNumIters"]

    dirpaths = mg_cfd_output_dirpaths

    run_id = 0

    agg_dfs = {k:None for k in cats}

    for dp in dirpaths:
        for root, dirnames, filenames in os.walk(dp):
            data_found = False
            for cat in cats:
                for filename in fnmatch.filter(filenames, cat+'.csv'):
                    data_found = True
            if data_found:
                run_id += 1
                for cat in cats:
                    for filename in fnmatch.filter(filenames, cat+'.csv'):
                        df_filepath = os.path.join(root, filename)
                        df = clean_pd_read_csv(df_filepath)

                        if cat == "Attributes":
                            simd_failed = did_simd_fail(root, infer_field(root, "CC")) and (infer_field(root, "SIMD")=="Y")
                            df = df.append({"Attribute":"SIMD failed", "Value":simd_failed}, ignore_index=True)

                        df['Run ID'] = run_id
                        # df = pd.pivot(df, index="Run ID", columns="Attribute", values="Value")

                        new_col_ordering_template = ["Run ID", "Attribute", "Event", "MG level", "Loop", "ThreadNum", "Instruction", "Value", "Count", "Time", "NumIters"]
                        new_col_ordering = [c for c in new_col_ordering_template if c in df.columns.values]
                        if len(new_col_ordering) != len(df.columns.values):
                            raise Exception("New column ordering is missing {0} columns: {1}", len(df.columns.values)-len(new_col_ordering), Set(df.columns.values).difference(new_col_ordering))
                        df = df[new_col_ordering]

                        if agg_dfs[cat] is None:
                            agg_dfs[cat] = df
                        else:
                            agg_dfs[cat] = safe_pd_append(agg_dfs[cat], df)
                        agg_dfs[cat] = agg_dfs[cat].sort_values(by=new_col_ordering)

    for cat in cats:
        if agg_dfs[cat] is None:
            print("NOTICE: Failed to find any '{0}.csv' output files to collate".format(cat))
            continue

    if "PAPI" in agg_dfs and not agg_dfs["PAPI"] is None:
        papi_df = agg_dfs["PAPI"]
        offcore_dram_read_events = [
            "OFFCORE_REQUESTS:ALL_DATA_READ",
            "OFFCORE_RESPONSE_[01]:(DMND_DATA_RD|PF_DATA_RD):(DMND_DATA_RD|PF_DATA_RD)",
            "OFFCORE_RESPONSE_[01]:(DMND_DATA_RD|PF_DATA_RD)"
            ]
        f = None
        for e in offcore_dram_read_events:
            g = papi_df["Event"].str.match(r''+e)
            if f is None:
                f = g
            else:
                f = np.logical_or(f, g)
        papi_df.loc[f,"Event"] = "dram read GB"
        papi_df.loc[f,"Count"] = papi_df.loc[f,"Count"] * 64 / 1e9
        agg_dfs["PAPI"] = papi_df

    # Drop job ID columns with just one value:
    att_df = agg_dfs["Attributes"]
    att_keys = att_df["Attribute"].drop_duplicates()
    invariant_atts = []
    variant_atts = []
    for k in att_keys:
        key_values = att_df.loc[att_df["Attribute"]==k, "Value"].drop_duplicates()
        if key_values.shape[0] == 1:
            invariant_atts.append(k)
        else:
            variant_atts.append(k)
    att_pruned_df = att_df.loc[att_df["Attribute"].isin(list(Set(variant_atts).union(essential_colnames)))]
    if att_pruned_df.shape[0] > 0:
        agg_dfs["Attributes"] = att_pruned_df

    att_df = agg_dfs["Attributes"]
    att_df = att_df.pivot(columns='Attribute', index='Run ID', values='Value')
    att_df = att_df.reset_index()

    invariant_atts = list(Set(invariant_atts).intersection(att_df.columns.values))
    variant_atts = list(Set(variant_atts).intersection(att_df.columns.values))

    ## Now combine attributes with the other tables:
    for cat in list(Set(cats).difference(["Attributes"])):
        if not agg_dfs[cat] is None:
            agg_dfs[cat] = agg_dfs[cat].merge(att_df)

            data_colnames = [c for c in agg_dfs[cat].columns.values if not c in att_df.columns.values]
            new_col_ordering = ["Run ID"] + invariant_atts + variant_atts + data_colnames
            if len(new_col_ordering) != len(agg_dfs[cat].columns.values):
                raise Exception("New column ordering has missing/additional columns")
            agg_dfs[cat] = agg_dfs[cat][new_col_ordering]

    for cat in cats:
        if not agg_dfs[cat] is None:
            if agg_dfs[cat].shape[0] == 0:
                raise Exception("Collated data table for '{0}' is empty".format(cat))
            agg_fp = os.path.join(prepared_output_dirpath, cat+".csv")
            if not os.path.isdir(prepared_output_dirpath):
                os.mkdir(prepared_output_dirpath)
            agg_dfs[cat].to_csv(agg_fp, index=False)

def aggregate():
    for cat in ["Times"]:
        df_filepath = os.path.join(prepared_output_dirpath,cat+".csv")
        if not os.path.isfile(df_filepath):
            continue
        print("Aggregating " + cat)
        df = clean_pd_read_csv(df_filepath)
        if "ThreadNum" in df.columns.values:
            df = df.drop("ThreadNum", axis=1)
        job_id_colnames = get_job_id_colnames(df)
        data_colnames = [c for c in df.columns.values if not c in job_id_colnames]
        if "Run ID" in job_id_colnames:
            job_id_colnames.remove("Run ID")

        df = df[job_id_colnames + data_colnames]
        df[data_colnames] = df[data_colnames].replace(0, np.NaN)

        df_agg = df.groupby(job_id_colnames)

        df_mean = df_agg.mean().reset_index().replace(np.NaN, 0.0)
        df_mean = df_mean[job_id_colnames + data_colnames]
        out_filepath = os.path.join(prepared_output_dirpath, cat+".mean.csv")
        df_mean.to_csv(out_filepath, index=False)

        if cat == "Times":
            # Calculate STDEV as % of mean:
            df_std = df_agg.std().reset_index().replace(np.NaN, 0.0)
            renames = {}
            for x in data_colnames:
                renames[x] = x+"_mean"
            df_std = df_std.merge(df_mean.rename(index=str, columns=renames))
            for c in data_colnames:
                df_std[c] = df_std[c] / df_std[c+"_mean"]
                df_std.loc[df_std[c+"_mean"]==0.0,c] = 0.0
                df_std.loc[df_std[c+"_mean"]==0,c]   = 0.0
            df_std = df_std.drop([x+"_mean" for x in data_colnames], axis=1)
            out_filepath = os.path.join(prepared_output_dirpath, cat+".std_pct.csv")
            df_std.to_csv(out_filepath, index=False)

    cats = ["instruction-counts", "instruction-counts-categorised"]
    for cat in cats:
        df_filepath = os.path.join(prepared_output_dirpath,cat+".csv")
        if os.path.isfile(df_filepath):
            print("Aggregating " + cat)
            df = clean_pd_read_csv(df_filepath)
            if "ThreadNum" in df.columns.values:
                df = df.drop("ThreadNum", axis=1)
            job_id_colnames = get_job_id_colnames(df)
            data_colnames = [c for c in df.columns.values if not c in job_id_colnames]
            if "Run ID" in job_id_colnames:
                job_id_colnames.remove("Run ID")

            df = df[job_id_colnames + data_colnames]
            df_agg = df.groupby(job_id_colnames)

            df_mean = df_agg.mean().reset_index()
            out_filepath = os.path.join(prepared_output_dirpath, cat+".csv")
            df_mean.to_csv(out_filepath, index=False)

    cat = "LoopNumIters"
    df_filepath = os.path.join(prepared_output_dirpath,cat+".csv")
    if os.path.isfile(df_filepath):
        print("Aggregating " + cat)
        df = clean_pd_read_csv(df_filepath)
        df["Metric"] = "#iterations"
        job_id_colnames = get_job_id_colnames(df)
        data_colnames = [c for c in df.columns.values if not c in job_id_colnames]
        if "Run ID" in job_id_colnames:
            job_id_colnames.remove("Run ID")

        df = df[job_id_colnames + data_colnames]
        df[data_colnames] = df[data_colnames].replace(0, np.NaN)

        df_agg = df.groupby(job_id_colnames)

        df_mean = df_agg.mean().reset_index().replace(np.NaN, 0.0)
        ## Next, compute sum and max across threads within each run:
        if "ThreadNum" in df.columns.values:
            del job_id_colnames[job_id_colnames.index("ThreadNum")]
            df_mean2 = df_mean.drop("ThreadNum", axis=1)
        else:
            df_mean2 = df_mean
        df_agg2 = df_mean2.groupby(job_id_colnames)
        df_sum = df_agg2.sum().reset_index()
        df_sum["Metric"] = "#iterations_SUM"
        df_max = df_agg2.max().reset_index()
        df_max["Metric"] = "#iterations_MAX"
        df_agg3 = df_sum.append(df_max)
        out_filepath = os.path.join(prepared_output_dirpath, cat+".csv")
        df_agg3.to_csv(out_filepath, index=False)

    cat = "PAPI"
    papi_df_filepath = os.path.join(prepared_output_dirpath,cat+".csv")
    if os.path.isfile(papi_df_filepath):
        print("Aggregating " + cat)
        df = clean_pd_read_csv(papi_df_filepath)
        job_id_colnames = get_job_id_colnames(df)
        data_colnames = [c for c in df.columns.values if not c in job_id_colnames]
        if "Run ID" in job_id_colnames:
            job_id_colnames.remove("Run ID")

        df = df[job_id_colnames + data_colnames]
        df[data_colnames] = df[data_colnames].replace(0.0, np.NaN)

        df_grps = df.groupby(job_id_colnames)

        ## First, compute per-thread average across repeat runs:
        df_thread_means = df_grps.mean().reset_index().replace(np.NaN, 0.0)

        ## Next, compute sum and max across threads within each run:
        if "ThreadNum" in df.columns.values:
            job_id_colnames.remove("ThreadNum")
            df_thread_means = df_thread_means.drop("ThreadNum", axis=1)

        df_grps = df_thread_means.groupby(job_id_colnames)
        df_run_sums = df_grps.sum().reset_index()
        df_run_maxs = df_grps.max().reset_index()

        df_grps = df_thread_means.replace(0, np.NaN).groupby(job_id_colnames)
        df_run_means = df_grps.mean().reset_index().replace(np.NaN, 0.0)

        # Calculate STDEV as % of mean:
        df_std = df_grps.std().reset_index().replace(np.NaN, 0.0)
        df_std_pct = safe_frame_divide(df_std, df_run_means)

        for pe in Set(df_run_sums["Event"]):
            df_run_sums.loc[df_run_sums["Event"]==pe, "Event"] = pe
        for pe in Set(df_run_maxs["Event"]):
            df_run_maxs.loc[df_run_maxs["Event"]==pe, "Event"] = pe+".THREADS_MAX"
        for pe in Set(df_run_means["Event"]):
            df_run_means.loc[df_run_means["Event"]==pe, "Event"] = pe+".THREADS_MEAN"
        df_agg = safe_pd_append(safe_pd_append(df_run_sums, df_run_maxs), df_run_means)

        out_filepath = os.path.join(prepared_output_dirpath, cat+".mean.csv")
        df_agg.to_csv(out_filepath, index=False)

        out_filepath = os.path.join(prepared_output_dirpath, cat+".std_pct.csv")
        df_std_pct.to_csv(out_filepath, index=False)

def count_fp_ins():
    insn_df_filepath = os.path.join(prepared_output_dirpath, "instruction-counts.csv")
    if not os.path.isfile(insn_df_filepath):
        return None

    intel_fp_insns = [
        "[v]?[u]?comisd", "[v]?add[sp][sd]", "[v]?maxsd", 
        "[v]?min[sp]d", "[v]?sub[sp][sd]", "vcmpltsd", 
        "[v]?div[sp][sd]", "[v]?sqrt[sp][sd]", 
        "vrcp14[sp]d", "vrcp28[sp]d", 
        "vrsqrt14[sp]d", "vrsqrt28[sp]d", 
        "[v]?mul[sp][sd]", "vgetexpsd", "vgetmantsd", "vscalefsd"
    ]
    intel_fp_fma_insns = ["vf[n]?madd[0-9]*[sp][sd]", "vf[n]?msub[0-9]*[sp][sd]"]

    aarch64_fp_insns = [
        "[v]?fadd", "[v]?fsub", 
        "[v]?fdiv", "[v]?fsqrt", 
        "[v]?fmul",
        "vfneg"
    ]
    aarch64_fp_fma_insns = [
        "f[n]?madd", "f[n]?msub", "vfml[as]"
    ]

    insn_df = clean_pd_read_csv(os.path.join(insn_df_filepath))
    fp_df = safe_pd_filter(insn_df, "Loop", key_loops)

    fp_df["FLOPs/iter"] = 0
    fp_df["FP ins/iter"] = 0
    for f in intel_fp_insns:
        f = fp_df["Instruction"].str.match(f)
        fp_df.loc[f,"FLOPs/iter"]  = fp_df.loc[f,"Count"]
        fp_df.loc[f,"FP ins/iter"] = fp_df.loc[f,"Count"]
    for f in intel_fp_fma_insns:
        f = fp_df["Instruction"].str.match(f)
        fp_df.loc[f,"FLOPs/iter"]  = fp_df.loc[f,"Count"]*2
        fp_df.loc[f,"FP ins/iter"] = fp_df.loc[f,"Count"]
    for f in aarch64_fp_insns:
        f = fp_df["Instruction"].str.match(f)
        fp_df.loc[f,"FLOPs/iter"]  = fp_df.loc[f,"Count"]
        fp_df.loc[f,"FP ins/iter"] = fp_df.loc[f,"Count"]
    for f in aarch64_fp_fma_insns:
        f = fp_df["Instruction"].str.match(f)
        fp_df.loc[f,"FLOPs/iter"]  = fp_df.loc[f,"Count"]*2
        fp_df.loc[f,"FP ins/iter"] = fp_df.loc[f,"Count"]

    fp_df = fp_df.drop(["Instruction", "Count"], axis=1)
    grp_colnames = [c for c in fp_df.columns.values if not c in ["FLOPs/iter", "FP ins/iter"] ]
    fp_df_grps = fp_df.groupby(grp_colnames)
    fp_df_sum = fp_df_grps.sum().reset_index().replace(np.NaN, 0.0)

    if "SIMD len" in fp_df_sum.columns.values:
        simd_mask = np.invert(fp_df_sum["SIMD failed"])
        fp_df_sum.loc[simd_mask, "FLOPs/iter"] *= fp_df_sum.loc[simd_mask, "SIMD len"]

    return fp_df_sum

def combine_all():
    data_all = None
    metrics = []

    times_df_filepath = os.path.join(prepared_output_dirpath, "Times.mean.csv")
    if os.path.isfile(times_df_filepath):
        times_df = clean_pd_read_csv(times_df_filepath)
        times_df["Metric"] = "Runtime"
        times_df = times_df.rename(index=str, columns={"Time":"Value"})
        data_colnames = get_data_colnames(times_df)
        data_all = times_df
        metrics.append("Runtime")

    papi_df_filepath = os.path.join(prepared_output_dirpath, "PAPI.mean.csv")
    if os.path.isfile(papi_df_filepath):
        papi_df = clean_pd_read_csv(papi_df_filepath)
        metrics += list(Set(papi_df["Event"]))
        papi_df = papi_df.rename(index=str, columns={"Event":"Metric", "Count":"Value"})

        if data_all is None:
            data_all = papi_df
        else:
            data_all = safe_pd_append(data_all, papi_df)

    if data_all is None:
        return

    iters_df_filepath = os.path.join(prepared_output_dirpath, "LoopNumIters.csv")
    if os.path.isfile(iters_df_filepath):
        iters_df = clean_pd_read_csv(iters_df_filepath)
        iters_df = safe_pd_filter(iters_df, "Metric", "#iterations_SUM")
        if "SIMD len" in iters_df.columns.values:
            ## Need the number of SIMD iterations. Currently 'iterations_SUM'
            ## counts the number of serial iterations.
            simd_mask = np.invert(iters_df["SIMD failed"])
            data_colnames = get_data_colnames(iters_df)
            for dc in data_colnames:
                iters_df.loc[simd_mask, dc] /= iters_df.loc[simd_mask, "SIMD len"]
        flux_iters_df = safe_pd_filter(iters_df, "Loop", "compute_flux_edge")

        fp_counts = count_fp_ins()
        if not fp_counts is None:
            ## Calculate GFLOPs
            flops_total = fp_counts.drop("FP ins/iter", axis=1).merge(flux_iters_df)
            flops_total["Metric"] = "GFLOPs"
            flops_total["Value"] = flops_total["FLOPs/iter"] * flops_total["NumIters"] / 1e9
            flops_total.drop(["NumIters", "FLOPs/iter"], axis=1, inplace=True)
            if data_all is None:
                data_all = flops_total
            else:
                data_all = safe_pd_append(data_all, flops_total)
            metrics.append("GFLOPs")

            if "PAPI_TOT_CYC.THREADS_MAX" in metrics:
                ## Calculate IPC of FP instructions:
                fp_total = fp_counts.drop("FLOPs/iter", axis=1).merge(flux_iters_df)
                fp_total["Metric"] = "FP total"
                fp_total["Value"] = fp_total["FP ins/iter"] * fp_total["NumIters"]
                fp_total.drop(["NumIters", "FP ins/iter"], axis=1, inplace=True)

                cyc_data = safe_pd_filter(data_all, "Metric", "PAPI_TOT_CYC")
                cyc_data = safe_pd_filter(cyc_data, "Loop", key_loops)
                fp_ipc = safe_frame_divide(fp_total, cyc_data)
                fp_ipc["Metric"] = "FP IPC"
                data_all = safe_pd_append(data_all, fp_ipc)
                metrics.append("FP IPC")

            ## Append FP ins/iter:
            fpInsPerIter = fp_counts.drop("FLOPs/iter", axis=1)
            fpInsPerIter["Metric"] = "FP ins/iter"
            fpInsPerIter = fpInsPerIter.rename(index=str, columns={"FP ins/iter":"Value"})
            # Metric 'FP ins/iter' is invariant to MG level, but appending to 'data_all' requires it, 
            # so need to add a level column:
            fpInsPerIter["MG level"] = 0
            data_all = safe_pd_append(data_all, fpInsPerIter)

    if "Runtime" in metrics and "dram read GB" in metrics:
        ## Calculate GB/sec
        runtime_data = data_all[data_all["Metric"]=="Runtime"]
        gb_data = data_all[data_all["Metric"]=="dram read GB"]
        gbsec_data = safe_frame_divide(gb_data, runtime_data.drop("Metric", axis=1))
        gbsec_data["Metric"] = "dram read GB/sec"
        data_all = safe_pd_append(data_all, gbsec_data)
    if "Runtime" in metrics and "GFLOPs" in metrics:
        ## Calculate GFLOPs/sec
        runtime_data = data_all[data_all["Metric"]=="Runtime"]
        gflops_data = data_all[data_all["Metric"]=="GFLOPs"]
        gflops_data = safe_frame_divide(gflops_data, runtime_data.drop("Metric", axis=1))
        gflops_data["Metric"] = "GFLOPs/sec"
        data_all = safe_pd_append(data_all, gflops_data)
    if "GFLOPs" in metrics and "dram read GB" in metrics:
        ## Calculate Flops/Byte
        gflops_data = data_all[data_all["Metric"]=="GFLOPs"]
        gb_data = data_all[data_all["Metric"]=="dram read GB"]
        arith_intensity_data = safe_frame_divide(gflops_data, gb_data.drop("Metric", axis=1))
        arith_intensity_data["Metric"] = "Flops/Byte"
        data_all = safe_pd_append(data_all, arith_intensity_data)
    if "Runtime" in metrics and "PAPI_TOT_CYC.THREADS_MAX" in metrics:
        ## Calculate GHz
        cyc_data = data_all[data_all["Metric"]=="PAPI_TOT_CYC.THREADS_MAX"]
        runtime_data = data_all[data_all["Metric"]=="Runtime"]
        ghz = safe_frame_divide(cyc_data, runtime_data.drop("Metric", axis=1))
        ghz["Value"] /= 1e9
        ghz["Metric"] = "GHz"
        data_all = safe_pd_append(data_all, ghz)
    if "GFLOPs" in metrics and "PAPI_TOT_CYC.THREADS_MAX" in metrics:
        ## Calculate flops/cycle
        # cyc_data = data_all[data_all["metric"]=="PAPI_TOT_CYC.THREADS_MAX"]
        cyc_data = data_all[data_all["Metric"]=="PAPI_TOT_CYC"]
        gflops_data = data_all[data_all["Metric"]=="GFLOPs"]
        flops_per_cyc = safe_frame_divide(gflops_data, cyc_data.drop("Metric", axis=1))
        flops_per_cyc["Value"] *= 1e9
        flops_per_cyc["Metric"] = "Flops/Cycle"
        data_all = safe_pd_append(data_all, flops_per_cyc)

    # Drop runtime:
    f = data_all["Metric"] == "Runtime"
    data_all = data_all[np.logical_not(f)]

    # Drop PAPI events:
    f = data_all["Metric"].str.contains("PAPI_")
    data_all = data_all[np.logical_not(f)]
    f = data_all["Metric"].str.contains("OFFCORE_")
    data_all = data_all[np.logical_not(f)]

    # Drop intermediate derivations from PAPI events:
    data_all = data_all[data_all["Metric"]!="GFLOPs"]
    data_all = data_all[data_all["Metric"]!="dram read GB"]
    data_all = data_all[data_all["Metric"]!="dram read GB.THREADS_MAX"]
    data_all = data_all[data_all["Metric"]!="dram read GB.THREADS_MEAN"]

    # Round to N significant figures:
    N=3
    # N=2
    f = data_all["Value"]!=0.0
    data_all.loc[f,"Value"] = data_all.loc[f,"Value"].apply(lambda x: round(x, N-1-int(floor(log10(abs(x))))))

    # Reorder columns like so: <job id columns>, <flux() columns>, <all other kernel columns>
    data_colnames = get_data_colnames(data_all)
    data_all = data_all[get_job_id_colnames(data_all)+data_colnames]

    if "Flux options" in data_all.columns.values:
        data_all = data_all[data_all["Flux options"]==""]
        data_all.drop("Flux options", axis=1, inplace=True)

    if "Flux variant" in data_all.columns.values:
        data_all = data_all[data_all["Flux variant"]=="Normal"]
        data_all.drop("Flux variant", axis=1, inplace=True)

    if not data_all is None:
        data_all.to_csv(os.path.join(prepared_output_dirpath, "Performance-Metrics-(derived).csv"), index=False)

if not assembly_analyser_dirpath is None:
    analyse_object_files()
collate_csvs()
aggregate()

combine_all()

print("Aggregated data written to folder: " + prepared_output_dirpath)
