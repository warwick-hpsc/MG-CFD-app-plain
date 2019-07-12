import os
import pandas as pd
import numpy as np
from sets import Set
import fnmatch
import argparse
import inspect

script_dirpath = os.path.dirname(os.path.realpath(__file__))
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

compile_info = {}
compile_info["compiler"] = "intel"
compile_info["SIMD len"] = 1
ins_per_iter = -1

kernels = ["flux", "update", "compute_step", "time_step", "up", "down", "indirect_rw"]

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
        if v.startswith("insn."):
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

            # print("Processing: " + output_dirpath)

            kernel_to_object = {}

            log_filename = ""
            log_filename_candidates = [x+".stdout" for x in ["sbatch", "moab", "lsf", "pbs"]]
            log_filename_candidates.append("submit.log")
            for lfc in log_filename_candidates:
                if os.path.isfile(os.path.join(output_dirpath, lfc)):
                    log_filename = lfc
                    break
            if log_filename == "":
                raise IOError("Cannot find a log file for run: " + output_dirpath)

            if grep("-DFLUX_CRIPPLE", os.path.join(output_dirpath, log_filename)):
                kernel_to_object["compute_flux_edge_crippled"] = "flux_loops.o"
            else:
                kernel_to_object["compute_flux_edge"] = "flux_loops.o"
            kernel_to_object["indirect_rw"] = "indirect_rw_loop.o"

            loops_tally_df = None
            for k in kernel_to_object.keys():
                obj_filepath = os.path.join(output_dirpath, "objects", kernel_to_object[k])
                loop, asm_loop_filepath = extract_loop_kernel_from_obj(obj_filepath, compile_info, ins_per_iter, k)
                loop_tally = count_loop_instructions(asm_loop_filepath, loop)
                for i in loop_tally.keys():
                    loop_tally["insn."+i] = loop_tally[i]
                    del loop_tally[i]

                if k == "compute_flux_edge_crippled":
                    k = "compute_flux_edge"
                loop_tally["kernel"] = k

                if loops_tally_df is None:
                    loops_tally_df = pd.DataFrame.from_dict({k:[v] for k,v in loop_tally.iteritems()})
                else:
                    for f in Set(loops_tally_df.keys()).difference(Set(loop_tally.keys())):
                        loop_tally[f] = 0
                    for f in Set(loop_tally.keys()).difference(Set(loops_tally_df.keys())):
                        loops_tally_df[f] = 0

                    if "sort" in inspect.getargspec(pd.DataFrame.append)[0]:
                        loops_tally_df = loops_tally_df.append(pd.DataFrame.from_dict({k:[v] for k,v in loop_tally.iteritems()}), sort=True)
                    else:
                        loops_tally_df = loops_tally_df.append(pd.DataFrame.from_dict({k:[v] for k,v in loop_tally.iteritems()}))

            job_id_df = get_output_run_config(output_dirpath)
            df = job_id_df.join(loops_tally_df)
            df.to_csv(ic_filepath, index=False)

def collate_csvs():
    cats = ["Times", "PAPI", "instruction-counts", "LoopNumIters"]

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
                            df_data_col_names = get_data_colnames(df_agg)
                            for d in df_missing_cols:
                                if d in df_data_col_names:
                                    df[d] = 0

                        df_agg_missing_cols = Set(df.columns.values).difference(Set(df_agg.columns.values))
                        if len(df_agg_missing_cols) > 0:
                            df_agg_data_col_names = get_data_colnames(df)
                            for d in df_agg_missing_cols:
                                if d in df_agg_data_col_names:
                                    df_agg[d] = 0

                        if "sort" in inspect.getargspec(pd.DataFrame.append)[0]:
                            df_agg = df_agg.append(df, sort=True)
                        else:
                            df_agg = df_agg.append(df)

        if df_agg is None:
            print("WARNING: Failed to find any '{0}' output files to collates".format(cat))
            continue

        if cat == "instruction-counts":
            df_agg = df_agg.drop("Size", axis=1)
        else:
            df_agg["Size_scale_factor"] = max(df_agg["Size"]) / df_agg["Size"]
            for dc in get_data_colnames(df_agg):
                df_agg[dc] *= df_agg["Size_scale_factor"]
            df_agg = df_agg.drop("Size", axis=1)
            df_agg = df_agg.drop("Size_scale_factor", axis=1)

        agg_fp = os.path.join(prepared_output_dirpath,cat+".csv")
        if not os.path.isdir(prepared_output_dirpath):
            os.mkdir(prepared_output_dirpath)
        df_agg.to_csv(agg_fp, index=False)

def aggregate():
    for cat in ["Times", "instruction-counts", "LoopNumIters"]:
        print("Aggregating " + cat)
        df_filepath = os.path.join(prepared_output_dirpath,cat+".csv")
        if not os.path.isfile(df_filepath):
            continue
        df = clean_pd_read_csv(df_filepath)
        if "ThreadNum" in df.columns.values:
            df = df.drop("ThreadNum", axis=1)
        job_id_colnames = get_job_id_colnames(df)
        df_agg = df.groupby(get_job_id_colnames(df), as_index=False)
        df_mean = df_agg.mean()
        out_filepath = os.path.join(prepared_output_dirpath, cat+".mean.csv")
        df_mean.to_csv(out_filepath, index=False)

    cat = "PAPI"
    print("Aggregating " + cat)
    df = clean_pd_read_csv(os.path.join(prepared_output_dirpath,cat+".csv"))
    job_id_colnames = get_job_id_colnames(df)
    ## First, compute per-thread average across repeat runs:
    df_agg = df.groupby(get_job_id_colnames(df), as_index=False)
    df_mean = df_agg.mean()
    ## Next, compute sum and max across threads within each run:
    del job_id_colnames[job_id_colnames.index("ThreadNum")]
    df_mean = df_mean.drop("ThreadNum", axis=1)
    df_agg2 = df_mean.groupby(job_id_colnames, as_index=False)
    df_sum = df_agg2.sum()
    df_max = df_agg2.max()
    for pe in Set(df_max["PAPI counter"]):
        df_max.loc[df_max["PAPI counter"]==pe, "PAPI counter"] = pe+"_MAX"
    df_agg3 = df_sum.append(df_max)
    out_filepath = os.path.join(prepared_output_dirpath, cat+".mean.csv")
    df_agg3.to_csv(out_filepath, index=False)

if not assembly_analyser_dirpath is None:
    analyse_object_files()
collate_csvs()
aggregate()

print("Aggregated data written to folder: " + prepared_output_dirpath)
