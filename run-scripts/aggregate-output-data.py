import os
import pandas as pd
import numpy as np
from sets import Set
import fnmatch
import argparse

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

kernels = ["flux", "update", "compute_step", "time_step", "up", "down", "indirect_rw"]

essential_colnames = ["CPU", "PAPI counter"]

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
        simd_possible = False
        if "FLUX_FISSION" in iters.loc[0,"Flux options"]:
            simd_possible = True
        elif "SIMD conflict avoidance strategy" in iters.columns.values and iters.loc[0,"SIMD conflict avoidance strategy"] != "" and iters.loc[0,"SIMD conflict avoidance strategy"] != "None":
            simd_possible = True
        if simd_possible:
            simd_len = iters.loc[0,"SIMD len"]
        else:
            simd_len = 1

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
            # 'num_iters' is counting non-vectorised iterations. If code 
            # was vectorised then fewer loop iterations will have been 
            # performed. Make that adjustment:
            num_iters /= simd_len

            # print("Kernel {0}, ins = {1}, iters = {2}, SIMD len = {3}".format(kernel, num_insn, num_iters, simd_len))

            ins_per_iter = float(num_insn) / float(num_iters)

    return ins_per_iter

def infer_compiler(output_dirpath):
    times_filepath = os.path.join(output_dirpath, "Times.csv")
    if not os.path.isfile(times_filepath):
        return "UNKNOWN"
    else:
        times = pd.read_csv(times_filepath)
        return times.loc[0, "CC"]

def infer_simd_len(output_dirpath):
    times_filepath = os.path.join(output_dirpath, "Times.csv")
    if not os.path.isfile(times_filepath):
        return 1
    else:
        times = pd.read_csv(times_filepath)
        return times.loc[0, "SIMD len"]

def analyse_object_files():
    print("Analysing object files")

    dirpaths = mg_cfd_output_dirpaths

    for dp in dirpaths:
        for o in os.listdir(dp):
            output_dirpath = os.path.join(dp, o)
            if not os.path.isdir(output_dirpath):
                continue

            ic_filepath = os.path.join(output_dirpath, "instruction-counts.csv")
            if os.path.isfile(ic_filepath):
                continue
            ic_cat_filepath = os.path.join(output_dirpath, "instruction-counts-categorised.csv")

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

            compile_info["compiler"] = infer_compiler(output_dirpath)
            compile_info["SIMD len"] = infer_simd_len(output_dirpath)

            loops_tally_df = None
            for k in kernel_to_object.keys():
                ins_per_iter = calc_ins_per_iter(output_dirpath, k)

                obj_filepath = os.path.join(output_dirpath, "objects", kernel_to_object[k])
                try:
                    loop, asm_loop_filepath = extract_loop_kernel_from_obj(obj_filepath, compile_info, ins_per_iter, k)
                except:
                    continue
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
                    ## The 'indirect rw' loop does not actually have any register spilling, so 
                    ## any that were identified by 'assembly-loop-extractor' are false positives. 
                    ## Assume that the same number of false positives are identified in the 'flux' kernel, 
                    ## so remove those:
                    claimed_indirect_rw_spills = loops_tally_df[loops_tally_df["kernel"]=="indirect_rw"].loc[0,"insn.LOAD_SPILLS"]
                    loops_tally_df["insn.LOAD_SPILLS"] -= claimed_indirect_rw_spills

                job_id_df = get_output_run_config(output_dirpath)
                df = job_id_df.join(loops_tally_df)
                df.to_csv(ic_filepath, index=False)

                f = categorise_aggregated_instructions_tally(ic_filepath)
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

                        df_agg = df_agg.append(df, sort=True)

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
        df_agg = df.groupby(get_job_id_colnames(df))

        df_mean = df_agg.mean().reset_index()
        out_filepath = os.path.join(prepared_output_dirpath, cat+".mean.csv")
        df_mean.to_csv(out_filepath, index=False)

        if cat == "Times":
            # Calculate STDEV as % of mean:
            df_std = df_agg.std().reset_index()
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
            df_agg = df.groupby(get_job_id_colnames(df))

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
        df_agg = df.groupby(get_job_id_colnames(df))

        df_mean = df_agg.mean().reset_index()
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
        ## First, compute per-thread average across repeat runs:
        df_agg = df.groupby(get_job_id_colnames(df))

        df_mean = df_agg.mean().reset_index()
        ## Next, compute sum and max across threads within each run:
        if "ThreadNum" in df.columns.values:
            del job_id_colnames[job_id_colnames.index("ThreadNum")]
            df_mean2 = df_mean.drop("ThreadNum", axis=1)
        else:
            df_mean2 = df_mean
        df_agg2 = df_mean2.groupby(job_id_colnames)
        df_sum = df_agg2.sum().reset_index()
        for pe in Set(df_sum["PAPI counter"]):
            df_sum.loc[df_sum["PAPI counter"]==pe, "PAPI counter"] = pe+"_SUM"
        df_max = df_agg2.max().reset_index()
        for pe in Set(df_max["PAPI counter"]):
            df_max.loc[df_max["PAPI counter"]==pe, "PAPI counter"] = pe+"_MAX"
        df_mean2 = df_agg2.mean().reset_index()
        for pe in Set(df_mean2["PAPI counter"]):
            df_mean2.loc[df_mean2["PAPI counter"]==pe, "PAPI counter"] = pe+"_MEAN"
        df_agg3 = df_sum.append(df_max, sort=True).append(df_mean2, sort=True)
        out_filepath = os.path.join(prepared_output_dirpath, cat+".mean.csv")
        df_agg3.to_csv(out_filepath, index=False)

        # Calculate STDEV as % of mean:
        df_std = df_agg.std().reset_index()
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

def count_flops():
    insn_df_filepath = os.path.join(prepared_output_dirpath, "instruction-counts.csv")
    if not os.path.isfile(insn_df_filepath):
        return 0

    intel_fp_insns = [
        "[v]?add[sp][sd]", "[v]?[u]?comisd", "[v]?maxsd", 
        "[v]?min[sp]d", "[v]?sub[sp][sd]", "[v]?div[sp][sd]", 
        "[v]?sqrt[sp][sd]", "vrsqrt14[sp]d", "vrsqrt28[sp]d", 
        "[v]?mul[sp][sd]", 
        "[v]?dpp[sd]", "[v]?movhp[sd]", "[v]?unpck[hl]p[sd]",
        "[v]?xorp[sd]", "vandpd", "vbroadcasts[sd]",
        "vextractf128", "vinsertf128", "vmovddup", "vpermilpd"]
    intel_fp_fma_insns = ["vf[n]?madd[0-9]*[sp][sd]", "vf[n]?msub[0-9]*[sp][sd]"]
    insn_df = clean_pd_read_csv(os.path.join(insn_df_filepath))
    insn_df = insn_df[insn_df["kernel"]=="compute_flux_edge"]
    if "Num threads" in insn_df.columns.values:
        insn_df = insn_df[insn_df["Num threads"] == insn_df["Num threads"].min()]
    insn_df["FLOPs"] = 0
    for colname in insn_df.columns.values:
        if colname.startswith("insn."):
            insn = colname.replace("insn.", "")
            fp_detected = False
            for f in intel_fp_insns:
                if re.match(f, insn):
                    insn_df["FLOPs"] += insn_df[colname]
                    fp_detected = True
                    break
            if not fp_detected:
                for f in intel_fp_fma_insns:
                    if re.match(f, insn):
                        insn_df["FLOPs"] += 2*insn_df[colname]
                        fp_detected = True
                        break

            insn_df = insn_df.drop(colname, axis=1)

    return insn_df.loc[0, "FLOPs"]

def combine_all():
    data_all = None

    times_df_filepath = os.path.join(prepared_output_dirpath, "Times.mean.csv")
    if os.path.isfile(times_df_filepath):
        times_df = clean_pd_read_csv(times_df_filepath)
        times_df["counter"] = "runtime"
        data_all = times_df

    papi_df_filepath = os.path.join(prepared_output_dirpath, "PAPI.mean.csv")
    if os.path.isfile(papi_df_filepath):
        papi_df = clean_pd_read_csv(papi_df_filepath)
        papi_df = papi_df[papi_df["PAPI counter"].str.contains("_SUM")]
        papi_df = papi_df.rename(index=str, columns={"PAPI counter":"counter"})

        if data_all is None:
            data_all = papi_df
        else:
            data_all = data_all.append(papi_df, sort=True)

    iters_df_filepath = os.path.join(prepared_output_dirpath, "LoopNumIters.csv")
    if os.path.isfile(iters_df_filepath):
        iters_df = clean_pd_read_csv(iters_df_filepath)
        iters_df = iters_df[iters_df["counter"].str.contains("_SUM")]

        flops_per_iter = count_flops()
        flops_total = iters_df[iters_df["counter"]=="#iterations_SUM"].copy()
        flops_total["counter"] = "GFLOPs"
        data_colnames = get_data_colnames(iters_df)
        flux_data_colnames = [c for c in data_colnames if c.startswith("flux")]
        flops_total[flux_data_colnames] = (flops_total[flux_data_colnames] / 1e9) * flops_per_iter
        other_data_colnames = [c for c in data_colnames if not c.startswith("flux")]
        flops_total[other_data_colnames] = 0.0
        if data_all is None:
            data_all = flops_total
        else:
            data_all = data_all.append(flops_total, sort=True)

    counters = list(Set(data_all["counter"]))
    if "runtime" in counters and "GB_SUM" in counters:
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
        gb_data["counter"] = "GB_SUM/sec"
        data_all = data_all.append(gb_data, sort=True)
    if "runtime" in counters and "GFLOPs" in counters:
        runtime_data = data_all[data_all["counter"]=="runtime"]
        gflops_data = data_all[data_all["counter"]=="GFLOPs"]
        data_colnames = get_data_colnames(runtime_data)
        renames = {}
        for x in data_colnames:
            renames[x] = x+"_runtime"
        runtime_data = runtime_data.rename(index=str, columns=renames)
        runtime_data = runtime_data.drop("counter", axis=1)
        gflops_data = gflops_data.drop("counter", axis=1)
        gflops_data = gflops_data.merge(runtime_data, validate="one_to_one")
        for cn in data_colnames:
            f = gflops_data[cn+"_runtime"] != 0.0
            gflops_data.loc[f,cn] /= gflops_data.loc[f,cn+"_runtime"]
            gflops_data = gflops_data.drop(cn+"_runtime", axis=1)
        gflops_data["counter"] = "GFLOPs/sec"
        data_all = data_all.append(gflops_data, sort=True)
    if "GFLOPs" in counters and "GB_SUM" in counters:
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
        gflops_data["counter"] = "Flops/Byte"
        data_all = data_all.append(gflops_data, sort=True)

    data_all = data_all[get_job_id_colnames(data_all)+get_data_colnames(data_all)]

    if not data_all is None:
        data_all.to_csv(os.path.join(prepared_output_dirpath, "all-data-combined.csv"), index=False)

if not assembly_analyser_dirpath is None:
    analyse_object_files()
collate_csvs()
aggregate()

combine_all()

print("Aggregated data written to folder: " + prepared_output_dirpath)
