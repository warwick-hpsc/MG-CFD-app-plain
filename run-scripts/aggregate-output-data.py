import os
import pandas as pd
import numpy as np
from sets import Set
import fnmatch
import argparse
import re
from math import floor, log10

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

kernels = ["flux", "update", "compute_step", "time_step", "restrict", "prolong", "indirect_rw"]

essential_colnames = ["CPU", "PAPI counter", "CC", "CC version", "Instruction set"]
essential_colnames += ["SIMD failed"]
essential_colnames += ["SIMD conflict avoidance strategy"]

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
    return top.append(bottom, sort=True)

def safe_pd_filter(df, field, value):
    if not field in df.columns.values:
        print("WARNING: field '{0}' not in df".format(field))
        return df
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
        else:
            is_kernel_col = False
            for k in kernels:
                for l in range(4):
                    if v == k+str(l):
                        is_kernel_col = True
                        break;
                if is_kernel_col:
                    break
            if is_kernel_col:
                data_colnames.append(v)
    return data_colnames

def get_job_id_colnames(df):
    job_id_colnames = list(Set(df.columns.values).difference(Set(get_data_colnames(df))))

    if len(job_id_colnames) == 0:
        print("ERROR: get_job_id_colnames() failed to find any.")
        sys.exit(-1)
    return job_id_colnames

def get_output_run_config(output_dirpath):
    times_df = clean_pd_read_csv(os.path.join(output_dirpath, "Times.csv"))

    job_id_columns = get_job_id_colnames(times_df)
    job_id_df = times_df[job_id_columns].iloc[0]
    job_id_df = job_id_df.to_frame()
    job_id_df = job_id_df.transpose()
    return job_id_df

def calc_ins_per_iter(output_dirpath, kernel):
    if kernel == "compute_flux_edge" or kernel == "compute_flux_edge_crippled":
        timer = "flux0"
    elif kernel == "indirect_rw":
        timer = "indirect_rw0"
    else:
        return -1

    ins_per_iter = -1

    papi_filepath = os.path.join(output_dirpath, "PAPI.csv")
    loop_num_iters_filepath = os.path.join(output_dirpath, "LoopNumIters.csv")

    if os.path.isfile(papi_filepath) and os.path.isfile(loop_num_iters_filepath):
        iters = clean_pd_read_csv(loop_num_iters_filepath)
        simd_failed = did_simd_fail(output_dirpath)
        if simd_failed:
            simd_len = 1
        else:
            simd_len = iters.loc[0,"SIMD len"]

        papi = clean_pd_read_csv(papi_filepath)
        if "PAPI_TOT_INS" in papi["PAPI counter"].unique():
            papi = papi[papi["PAPI counter"]=="PAPI_TOT_INS"]
            papi = papi.drop("PAPI counter", axis=1)

            job_id_colnames = get_job_id_colnames(papi)
            data_colnames = list(Set(papi.columns.values).difference(job_id_colnames))

            renames = {}
            for x in data_colnames:
                renames[x] = x+"_iters"
            iters = iters.rename(index=str, columns=renames)
            papi_and_iters = papi.merge(iters)
            if papi_and_iters.shape[0] != papi.shape[0]:
                raise Exception("merge of papi and iters failed")
            num_insn = papi_and_iters.loc[0, timer]

            num_iters = papi_and_iters.loc[0, timer+"_iters"]
            if num_iters == 0:
                ins_per_iter = 0.0
            else:
                num_iters = float(num_iters)
                # 'num_iters' is counting non-vectorised iterations. If code 
                # was vectorised then fewer loop iterations will have been 
                # performed. Make that adjustment:
                num_iters /= float(simd_len)

                # print("Kernel {0}, ins = {1}, iters = {2}, SIMD len = {3}".format(kernel, num_insn, num_iters, simd_len))

                ins_per_iter = float(num_insn) / float(num_iters)

    return ins_per_iter

def infer_field(output_dirpath, field):
    value = None

    times_filepath = os.path.join(output_dirpath, "Times.csv")
    if os.path.isfile(times_filepath):
        times = pd.read_csv(times_filepath)
        if field in times.columns.values:
            return times.loc[0, field]

    if field == "Precise FP":
        ## Infer this specific field from compiler stdout:
        log_filepath = os.path.join(output_dirpath, "submit.log")
        if not os.path.isfile(log_filepath):
            ## Must use old logic for deducing how this was compiled:
            ##   FP was precise for all runs except for AVX512 compiled with Intel, as 
            ##   Intel would segfault.
            if infer_field(output_dirpath, "CC")=="intel" and "AVX512" in infer_field(output_dirpath, "Instruction set"):
                value = "N"
            else:
                value = "Y"
        if grep("precise-fp", log_filepath):
            value = "Y"
        else:
            value = "N"

    elif field == "CC":
        ## Infer this specific field from make stdout:
        log_filepath = os.path.join(output_dirpath, "submit.log")
        if not os.path.isfile(log_filepath):
            raise Exception("Cannot infer field '{0}'".format(field))
        if grep("COMPILER=", log_filepath):
            if grep("COMPILER=clang", log_filepath):
                value = "clang"
            elif grep("COMPILER=gnu", log_filepath):
                value = "gnu"
            elif grep("COMPILER=intel", log_filepath):
                value = "intel"

    elif field == "SIMD":
        ## Infer this specific field from make stdout:
        log_filepath = os.path.join(output_dirpath, "submit.log")
        if not os.path.isfile(log_filepath):
            raise Exception("Cannot infer field '{0}'".format(field))
        if grep("-DSIMD-", log_filepath):
            value = True
        else:
            value = False

    elif field == "SIMD len":
        ## Infer this specific field from make stdout:
        log_filepath = os.path.join(output_dirpath, "submit.log")
        if not os.path.isfile(log_filepath):
            raise Exception("Cannot infer field '{0}'".format(field))
        if not grep("-DSIMD-", log_filepath):
            value = 1
        else:
            if grep("DBLS_PER_SIMD=2-", log_filepath):
                value = 2
            elif grep("DBLS_PER_SIMD=4-", log_filepath):
                value = 4
            elif grep("DBLS_PER_SIMD=8-", log_filepath):
                value = 8
            elif grep("DBLS_PER_SIMD=1-", log_filepath):
                value = 1

    elif field == "SIMD conflict avoidance strategy":
        ## Infer this specific field from make stdout:
        log_filepath = os.path.join(output_dirpath, "submit.log")
        if not os.path.isfile(log_filepath):
            raise Exception("Cannot infer field '{0}'".format(field))
        if grep("-DMANUAL_GATHER-", log_filepath) and grep("-DMANUAL_SCATTER-", log_filepath):
            value = "ManualGatherScatter"
        elif grep("-DMANUAL_GATHER-", log_filepath):
            value = "ManualGather"
        elif grep("-DMANUAL_SCATTER-", log_filepath):
            value = "ManualScatter"
        elif grep("-DCOLOURED_CONFLICT_AVOIDANCE-", log_filepath):
            if grep("-DBIN_COLOURED_VECTORS-", log_filepath):
                value = "ColouredEdgeVectors"
            elif grep("-DBIN_COLOURED_CONTIGUOUS-", log_filepath):
                value = "ColouredEdgesContiguous"
            else:
                value = "None"
        else:
            value = "None"

    else:
        times_filepath = os.path.join(output_dirpath, "Times.csv")
        if os.path.isfile(times_filepath):
            times = pd.read_csv(times_filepath)
            if not field in times.columns.values:
                raise Exception("Field '{0}' not in: {1}".format(field, output_dirpath))
            value = times.loc[0, field]

    if value is None:
        raise Exception("Cannot infer field '{0}'".format(field))
    
    return value

def did_simd_fail(output_dirpath, compiler=None):
    failed = False
    log_filepath =  os.path.join(output_dirpath, "flux_loops.o.log")

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

            kernel_to_object = {}

            run_filename_candidates = [x+".batch" for x in ["slurm", "pbs", "moab", "lsf"]]
            run_filename_candidates.append("run.sh")
            for rhc in run_filename_candidates:
                if os.path.isfile(os.path.join(output_dirpath, rhc)):
                    run_filename = rhc
                    break
            if run_filename == "":
                raise IOError("Cannot find a run script for run: " + output_dirpath)
            if grep("-DFLUX_CRIPPLE", os.path.join(output_dirpath, run_filename)):
                kernel_to_object["compute_flux_edge_crippled"] = "flux_loops.o"
            else:
                kernel_to_object["compute_flux_edge"] = "flux_loops.o"
            kernel_to_object["indirect_rw"] = "indirect_rw_loop.o"

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
            if "manual" in compile_info["SIMD CA scheme"].lower():
                if compile_info["compiler"] == "clang":
                    compile_info["manual CA block width"] = (compile_info["SIMD len"]*3)+2
                else:
                    compile_info["manual CA block width"] = compile_info["SIMD len"]

            loops_tally_df = None
            for k in kernel_to_object.keys():
                ins_per_iter = calc_ins_per_iter(output_dirpath, k)

                obj_filepath = os.path.join(output_dirpath, "objects", kernel_to_object[k])
                try:
                    asm_loop_filepath = extract_loop_kernel_from_obj(obj_filepath, compile_info, ins_per_iter, k, avx512cd_required, 10)
                except:
                    raise
                    # continue
                loop_tally = count_loop_instructions(asm_loop_filepath)
                for i in loop_tally.keys():
                    loop_tally["insn."+i] = loop_tally[i]
                    del loop_tally[i]

                if k == "compute_flux_edge_crippled":
                    k = "compute_flux_edge"
                loop_tally["kernel"] = k

                if loops_tally_df is None:
                    tmp_dict = {}
                    for k,v in loop_tally.iteritems():
                        tmp_dict[k] = [v]
                    loops_tally_df = pd.DataFrame.from_dict(tmp_dict)
                else:
                    for f in Set(loops_tally_df.keys()).difference(Set(loop_tally.keys())):
                        loop_tally[f] = 0
                    for f in Set(loop_tally.keys()).difference(Set(loops_tally_df.keys())):
                        loops_tally_df[f] = 0

                    tmp_dict = {}
                    for k,v in loop_tally.iteritems():
                        tmp_dict[k] = [v]
                    loops_tally_df = loops_tally_df.append(pd.DataFrame.from_dict(tmp_dict), sort=True)

            if not loops_tally_df is None:
                if "insn.LOAD_SPILLS" in loops_tally_df.columns.values:
                    if "indirect_rw" in loops_tally_df["kernel"]:
                        ## The 'indirect rw' loop does not actually have any register spilling, so 
                        ## any that were identified by 'assembly-loop-extractor' are false positives. 
                        ## Assume that the same number of false positives are identified in the 'flux' kernel, 
                        ## so remove those:
                        claimed_indirect_rw_spills = loops_tally_df[loops_tally_df["kernel"]=="indirect_rw"].loc[0,"insn.LOAD_SPILLS"]
                        loops_tally_df["insn.LOAD_SPILLS"] -= claimed_indirect_rw_spills

                job_id_df = get_output_run_config(output_dirpath)
                df = job_id_df.join(loops_tally_df)
                df.to_csv(ic_filepath, index=False)

                target_is_aarch64 = False
                raw_asm_filepath = os.path.join(output_dirpath, "objects", "flux_loops.o.raw-asm")
                for line in open(raw_asm_filepath):
                  if "aarch64" in line:
                    target_is_aarch64 = True
                    break

                if target_is_aarch64:
                    f = categorise_aggregated_instructions_tally_csv(ic_filepath, is_aarch64=True)
                else:
                    f = categorise_aggregated_instructions_tally_csv(ic_filepath, is_intel64=True)
                f.to_csv(ic_cat_filepath, index=False)

def collate_csvs():
    cats = ["Times", "PAPI", "instruction-counts", "instruction-counts-categorised", "LoopNumIters"]

    dirpaths = mg_cfd_output_dirpaths

    for cat in cats:
        print("Collating " + cat)
        df_agg = None
        for dp in dirpaths:
            for root, dirnames, filenames in os.walk(dp):
                for filename in fnmatch.filter(filenames, cat+'.csv'):
                    df_filepath = os.path.join(root, filename)
                    df = clean_pd_read_csv(df_filepath)

                    simd_failed = did_simd_fail(root, infer_field(root, "CC")) and (infer_field(root, "SIMD")=="Y")
                    df["SIMD failed"] = simd_failed

                    # Now require 'Precise FP' column:
                    if not "Precise FP" in df.columns.values:
                        df["Precise FP"] = infer_field(root, "Precise FP")

                    if df_agg is None:
                        df_agg = df
                    else:
                        df_missing_cols = Set(df_agg.columns.values).difference(Set(df.columns.values))
                        if len(df_missing_cols) > 0:
                            df_agg_data_col_names = get_data_colnames(df_agg)
                            for d in df_missing_cols:
                                if d in df_agg_data_col_names:
                                    df[d] = 0

                        df_agg_missing_cols = Set(df.columns.values).difference(Set(df_agg.columns.values))
                        if len(df_agg_missing_cols) > 0:
                            df_data_col_names = get_data_colnames(df)
                            for d in df_agg_missing_cols:
                                if d in df_data_col_names:
                                    df_agg[d] = 0

                        df_agg = safe_pd_append(df_agg, df)

        if df_agg is None:
            print("WARNING: Failed to find any '{0}' output files to collates".format(cat))
            continue

        if cat in ["instruction-counts", "instruction-counts-categorised"]:
            df_agg = df_agg.drop("Size", axis=1)
        else:
            df_agg["Size_scale_factor"] = max(df_agg["Size"]) / df_agg["Size"]
            for dc in get_data_colnames(df_agg):
                df_agg[dc] *= df_agg["Size_scale_factor"]
            df_agg = df_agg.drop("Size", axis=1)
            df_agg = df_agg.drop("Size_scale_factor", axis=1)

        if cat == "PAPI":
            f = df_agg["PAPI counter"]=="OFFCORE_RESPONSE_0:ANY_DATA:ANY_RESPONSE"
            data_colnames = get_data_colnames(df_agg)
            df_agg.loc[f,"PAPI counter"] = "GB"
            df_agg.loc[f,data_colnames] = df_agg.loc[f,data_colnames] * 64 / 1e9

        # Drop job ID columns with just one value:
        n_uniq = df_agg.apply(pd.Series.nunique)
        uniq_colnames = n_uniq[n_uniq==1].index
        job_id_columns = get_job_id_colnames(df_agg)
        uniq_colnames = list(Set(uniq_colnames).intersection(job_id_columns))
        uniq_colnames = list(Set(uniq_colnames).difference(essential_colnames))
        df_agg_clean = df_agg.drop(uniq_colnames, axis=1)
        if df_agg_clean.shape[1] > 0:
            df_agg = df_agg_clean

        agg_fp = os.path.join(prepared_output_dirpath,cat+".csv")
        if not os.path.isdir(prepared_output_dirpath):
            os.mkdir(prepared_output_dirpath)
        df_agg.to_csv(agg_fp, index=False)

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
        data_colnames = list(Set(df.columns.values).difference(job_id_colnames))

        job_id_colnames = sorted(job_id_colnames)
        data_colnames = sorted(data_colnames)
        df = df[job_id_colnames + data_colnames]

        df[data_colnames] = df[data_colnames].replace(0, np.NaN)

        df_agg = df.groupby(get_job_id_colnames(df))

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
            data_colnames = list(Set(df.columns.values).difference(job_id_colnames))

            job_id_colnames = sorted(job_id_colnames)
            data_colnames = sorted(data_colnames)
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
        df["counter"] = "#iterations"
        job_id_colnames = get_job_id_colnames(df)
        data_colnames = list(Set(df.columns.values).difference(job_id_colnames))

        df[data_colnames] = df[data_colnames].replace(0, np.NaN)

        job_id_colnames = sorted(job_id_colnames)
        data_colnames = sorted(data_colnames)
        df = df[job_id_colnames + data_colnames]

        df_agg = df.groupby(get_job_id_colnames(df))

        df_mean = df_agg.mean().reset_index().replace(np.NaN, 0.0)
        ## Next, compute sum and max across threads within each run:
        if "ThreadNum" in df.columns.values:
            del job_id_colnames[job_id_colnames.index("ThreadNum")]
            df_mean2 = df_mean.drop("ThreadNum", axis=1)
        else:
            df_mean2 = df_mean
        df_agg2 = df_mean2.groupby(job_id_colnames)
        df_sum = df_agg2.sum().reset_index()
        df_sum["counter"] = "#iterations_SUM"
        df_max = df_agg2.max().reset_index()
        df_max["counter"] = "#iterations_MAX"
        df_agg3 = df_sum.append(df_max)
        out_filepath = os.path.join(prepared_output_dirpath, cat+".csv")
        df_agg3.to_csv(out_filepath, index=False)

    cat = "PAPI"
    papi_df_filepath = os.path.join(prepared_output_dirpath,cat+".csv")
    if os.path.isfile(papi_df_filepath):
        print("Aggregating " + cat)
        df = clean_pd_read_csv(papi_df_filepath)
        job_id_colnames = get_job_id_colnames(df)
        data_colnames = list(Set(df.columns.values).difference(job_id_colnames))

        ## Exclude zero values from statistics:
        df[data_colnames] = df[data_colnames].replace(0.0, np.NaN)

        df_grps = df.groupby(get_job_id_colnames(df))

        ## First, compute per-thread average across repeat runs:
        df_thread_means = df_grps.mean().reset_index().replace(np.NaN, 0.0)

        ## Next, compute sum and max across threads within each run:
        if "ThreadNum" in df.columns.values:
            del job_id_colnames[job_id_colnames.index("ThreadNum")]
            df_thread_means = df_thread_means.drop("ThreadNum", axis=1)

        df_grps = df_thread_means.groupby(job_id_colnames)
        df_run_sums = df_grps.sum().reset_index()
        df_run_maxs = df_grps.max().reset_index()

        df_grps = df_thread_means.replace(0, np.NaN).groupby(job_id_colnames)
        df_run_means = df_grps.mean().reset_index().replace(np.NaN, 0.0)

        # Calculate STDEV as % of mean:
        df_std = df_grps.std().reset_index().replace(np.NaN, 0.0)
        df_std_pct = safe_frame_divide(df_std, df_run_means)

        for pe in Set(df_run_sums["PAPI counter"]):
            df_run_sums.loc[df_run_sums["PAPI counter"]==pe, "PAPI counter"] = pe+"_SUM"
        for pe in Set(df_run_maxs["PAPI counter"]):
            df_run_maxs.loc[df_run_maxs["PAPI counter"]==pe, "PAPI counter"] = pe+"_MAX"
        for pe in Set(df_run_means["PAPI counter"]):
            df_run_means.loc[df_run_means["PAPI counter"]==pe, "PAPI counter"] = pe+"_MEAN"
        df_agg = df_run_sums.append(df_run_maxs, sort=True).append(df_run_means, sort=True)
        
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
        "[v]?mul[sp][sd]"
    ]
    intel_fp_fma_insns = ["vf[n]?madd[0-9]*[sp][sd]", "vf[n]?msub[0-9]*[sp][sd]"]
    insn_df = clean_pd_read_csv(os.path.join(insn_df_filepath))
    insn_df = safe_pd_filter(insn_df, "kernel", "compute_flux_edge")
    insn_df["FLOPs/iter"] = 0
    insn_df["FP ins/iter"] = 0
    for colname in insn_df.columns.values:
        if colname.startswith("insn."):
            insn = colname.replace("insn.", "")
            fp_detected = False
            for f in intel_fp_insns:
                if re.match(f, insn):
                    insn_df["FLOPs/iter"]  += insn_df[colname]
                    insn_df["FP ins/iter"] += insn_df[colname]
                    fp_detected = True
                    break
            if not fp_detected:
                for f in intel_fp_fma_insns:
                    if re.match(f, insn):
                        insn_df["FLOPs/iter"]  += 2*insn_df[colname]
                        insn_df["FP ins/iter"] +=   insn_df[colname]
                        fp_detected = True
                        break

            insn_df = insn_df.drop(colname, axis=1)

    if "SIMD len" in insn_df.columns.values:
        simd_mask = np.invert(insn_df["SIMD failed"])
        insn_df.loc[simd_mask, "FLOPs/iter"] *= insn_df.loc[simd_mask, "SIMD len"]

    return insn_df

def combine_all():
    data_all = None

    times_df_filepath = os.path.join(prepared_output_dirpath, "Times.mean.csv")
    if os.path.isfile(times_df_filepath):
        times_df = clean_pd_read_csv(times_df_filepath)
        times_df["counter"] = "runtime"
        data_colnames = get_data_colnames(times_df)
        data_all = times_df

    papi_df_filepath = os.path.join(prepared_output_dirpath, "PAPI.mean.csv")
    if os.path.isfile(papi_df_filepath):
        papi_df = clean_pd_read_csv(papi_df_filepath)
        papi_df = papi_df.rename(index=str, columns={"PAPI counter":"counter"})

        if data_all is None:
            data_all = papi_df
        else:
            data_all = data_all.append(papi_df, sort=True)

    if data_all is None:
        return

    counters = list(Set(data_all["counter"]))
    data_colnames = get_data_colnames(data_all)
    flux_data_colnames = [c for c in data_colnames if c.startswith("flux")]
    other_data_colnames = [c for c in data_colnames if not c.startswith("flux")]

    iters_df_filepath = os.path.join(prepared_output_dirpath, "LoopNumIters.csv")
    if os.path.isfile(iters_df_filepath):
        iters_df = clean_pd_read_csv(iters_df_filepath)
        iters_df = safe_pd_filter(iters_df, "counter", "#iterations_SUM")
        if "SIMD len" in iters_df.columns.values:
            ## Need the number of SIMD iterations. Currently 'iterations_SUM'
            ## counts the number of serial iterations.
            simd_mask = np.invert(iters_df["SIMD failed"])
            for dc in data_colnames:
                iters_df.loc[simd_mask, dc] /= iters_df.loc[simd_mask, "SIMD len"]

        fp_counts = count_fp_ins()
        if not fp_counts is None:
            ## Calculate GFLOPs
            flops_total = safe_pd_merge(fp_counts.drop("FP ins/iter", axis=1), iters_df, "many_to_one")
            flops_total["counter"] = "GFLOPs"
            for f in flux_data_colnames:
                flops_total[f] = (flops_total[f] / 1e9) * flops_total["FLOPs/iter"]
            flops_total.drop("FLOPs/iter", axis=1, inplace=True)
            flops_total[other_data_colnames] = 0.0
            flops_total[data_colnames] = flops_total[data_colnames]
            if data_all is None:
                data_all = flops_total
            else:
                data_all = data_all.append(flops_total, sort=True)
            counters.append("GFLOPs")

            if "PAPI_TOT_CYC_MAX" in counters:
                ## Calculate IPC of FP instructions:
                fp_total = safe_pd_merge(fp_counts.drop("FLOPs/iter", axis=1), iters_df, "many_to_one")
                fp_total["counter"] = "FP total"
                flux_data_colnames = [c for c in data_colnames if c.startswith("flux")]
                for f in flux_data_colnames:
                    fp_total[f] = fp_total[f] * fp_total["FP ins/iter"]
                fp_total.drop("FP ins/iter", axis=1, inplace=True)
                fp_total[other_data_colnames] = 0.0
                fp_total[data_colnames] = fp_total[data_colnames]

                cyc_data = data_all[data_all["counter"]=="PAPI_TOT_CYC_SUM"]
                fp_ipc = safe_frame_divide(fp_total, cyc_data.drop("counter", axis=1))
                fp_ipc[data_colnames] = fp_ipc[data_colnames]
                fp_ipc["counter"] = "FP IPC"
                data_all = data_all.append(fp_ipc, sort=True)
                counters.append("FP IPC")

            ## Append FP ins/iter:
            fp_counts.drop("FLOPs/iter", axis=1, inplace=True)
            fp_counts["counter"] = "FP ins/iter"
            for dc in other_data_colnames:
                fp_counts[dc] = 0
            for dc in flux_data_colnames:
                fp_counts[dc] = fp_counts["FP ins/iter"].copy()
            fp_counts.drop("FP ins/iter", axis=1, inplace=True)
            data_all = safe_pd_append(data_all, fp_counts)

    if "runtime" in counters and "GB_SUM" in counters:
        ## Calculate GB_SUM/sec
        runtime_data = data_all[data_all["counter"]=="runtime"]
        gb_data = data_all[data_all["counter"]=="GB_SUM"]
        data_colnames = get_data_colnames(runtime_data)
        renames = {}
        for x in data_colnames:
            renames[x] = x+"_runtime"
        runtime_data = runtime_data.rename(index=str, columns=renames)
        runtime_data = runtime_data.drop("counter", axis=1)
        gb_data = gb_data.drop("counter", axis=1)
        gb_data = gb_data.merge(runtime_data, validate="one_to_one")
        for cn in data_colnames:
            f = gb_data[cn+"_runtime"] != 0.0
            gb_data.loc[f,cn] /= gb_data.loc[f,cn+"_runtime"]
            gb_data = gb_data.drop(cn+"_runtime", axis=1)
        ## Todo: test: the line below should produce the same result as block above
        # gb_data = safe_frame_divide(gb_data, runtime_data.drop("counter", axis=1))
        gb_data["counter"] = "GB_SUM/sec"
        data_all = data_all.append(gb_data, sort=True)
    if "runtime" in counters and "GFLOPs" in counters:
        ## Calculate GFLOPs/sec
        runtime_data = data_all[data_all["counter"]=="runtime"]
        gflops_data = data_all[data_all["counter"]=="GFLOPs"]
        gflops_data = safe_frame_divide(gflops_data, runtime_data.drop("counter", axis=1))
        gflops_data[data_colnames] = gflops_data[data_colnames]
        gflops_data["counter"] = "GFLOPs/sec"
        data_all = data_all.append(gflops_data, sort=True)
    if "GFLOPs" in counters and "GB_SUM" in counters:
        ## Calculate Flops/Byte
        gflops_data = data_all[data_all["counter"]=="GFLOPs"]
        gb_data = data_all[data_all["counter"]=="GB_SUM"]
        data_colnames = get_data_colnames(gflops_data)
        renames = {}
        for x in data_colnames:
            renames[x] = x+"_gb"
        gb_data = gb_data.rename(index=str, columns=renames)
        gb_data = gb_data.drop("counter", axis=1)
        gflops_data = gflops_data.drop("counter", axis=1)
        gflops_data = gflops_data.merge(gb_data, validate="one_to_one")
        for cn in data_colnames:
            f = gflops_data[cn+"_gb"] != 0.0
            gflops_data.loc[f,cn] /= gflops_data.loc[f,cn+"_gb"]
            f = gflops_data[cn+"_gb"] == 0.0
            gflops_data.loc[f,cn] = 0.0
            gflops_data = gflops_data.drop(cn+"_gb", axis=1)
        ## Todo: test: the line below should produce the same result as block above
        # gflops_data = safe_frame_divide(gflops_data, gb_data.drop("counter", axis=1))
        gflops_data["counter"] = "Flops/Byte"
        data_all = data_all.append(gflops_data, sort=True)
    if "runtime" in counters and "PAPI_TOT_CYC_MAX" in counters:
        ## Calculate GHz
        cyc_data = data_all[data_all["counter"]=="PAPI_TOT_CYC_MAX"]
        runtime_data = data_all[data_all["counter"]=="runtime"].copy()
        ghz = safe_frame_divide(cyc_data, runtime_data.drop("counter", axis=1))
        ghz[get_data_colnames(runtime_data)] /= 1e9
        ghz["counter"] = "GHz"
        ghz[data_colnames] = ghz[data_colnames]
        data_all = data_all.append(ghz, sort=True)
    if "GFLOPs" in counters and "PAPI_TOT_CYC_MAX" in counters:
        ## Calculate flops/cycle
        # cyc_data = data_all[data_all["counter"]=="PAPI_TOT_CYC_MAX"]
        cyc_data = data_all[data_all["counter"]=="PAPI_TOT_CYC_SUM"]
        gflops_data = data_all[data_all["counter"]=="GFLOPs"].copy()
        flops_per_cyc = safe_frame_divide(gflops_data, cyc_data.drop("counter", axis=1))
        flops_per_cyc[get_data_colnames(flops_per_cyc)] *= 1e9
        flops_per_cyc["counter"] = "Flops/Cycle"
        flops_per_cyc[data_colnames] = flops_per_cyc[data_colnames]
        data_all = data_all.append(flops_per_cyc, sort=True)

    # Drop PAPI counters:
    f = data_all["counter"].str.contains("PAPI_")
    data_all = data_all[np.logical_not(f)]

    # Drop GFLOPs:
    data_all = data_all[data_all["counter"]!="GFLOPs"]

    # Round to N significant figures:
    N=3
    # N=2
    data_colnames = get_data_colnames(data_all)
    for dc in data_colnames:
        f = data_all[dc]!=0.0
        data_all.loc[f,dc] = data_all.loc[f,dc].apply(lambda x: round(x, N-1-int(floor(log10(abs(x))))))

    # Reorder columns like so: <job id columns>, <flux() columns>, <all other kernel columns>
    data_colnames = get_data_colnames(data_all)
    flux_data_colnames = [c for c in data_colnames if c.startswith("flux")]
    other_data_colnames = [c for c in data_colnames if not c.startswith("flux")]
    data_all = data_all[get_job_id_colnames(data_all)+flux_data_colnames+other_data_colnames]

    if "Flux options" in data_all.columns.values:
        data_all = data_all[data_all["Flux options"]==""]
        data_all.drop("Flux options", axis=1, inplace=True)

    if "Flux variant" in data_all.columns.values:
        data_all = data_all[data_all["Flux variant"]=="Normal"]
        data_all.drop("Flux variant", axis=1, inplace=True)

    if not data_all is None:
        data_all.to_csv(os.path.join(prepared_output_dirpath, "all-data-combined.csv"), index=False)

if not assembly_analyser_dirpath is None:
    analyse_object_files()
collate_csvs()
aggregate()

combine_all()

print("Aggregated data written to folder: " + prepared_output_dirpath)
