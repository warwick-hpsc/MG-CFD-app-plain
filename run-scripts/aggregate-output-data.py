import os, shutil
import pandas as pd
import numpy as np
import fnmatch
import argparse
import re
import ast
from math import floor, log10, isnan

import sys
pyv = sys.version_info[0]
if pyv == 2:
    from sets import Set
elif pyv == 3:
    # Create alias:
    Set = set

# verbose = False
verbose = True

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
if os.path.isdir(prepared_output_dirpath):
    shutil.rmtree(prepared_output_dirpath)

if not assembly_analyser_dirpath is None:
    import imp
    imp.load_source('assembly_analysis', os.path.join(assembly_analyser_dirpath, "assembly_analysis.py"))
    from assembly_analysis import *
    imp.load_source('utils', os.path.join(assembly_analyser_dirpath, "utils.py"))
    from utils import *

compile_info = {}

essential_colnames = []
# essential_colnames += ["CPU", "CC", "CC version", "Instruction set"]
essential_colnames += ["CPU"] ## Performance model will need to know CPU
essential_colnames += ["Instruction set"] ## Performance model will need to know ISA
# essential_colnames += ["Event"]
# essential_colnames += ["SIMD failed", "SIMD conflict avoidance strategy"]
essential_colnames += ["SIMD failed"]
essential_colnames += ["ISA group"]

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
    simd_failed = do_logs_report_simd_fail(output_dirpath, kernel)
    ins_per_iter = -1

    papi_filepath = os.path.join(output_dirpath, "PAPI.csv")
    loop_num_iters_filepath = os.path.join(output_dirpath, "LoopNumIters.csv")

    if os.path.isfile(papi_filepath) and os.path.isfile(loop_num_iters_filepath):
        iters = clean_pd_read_csv(loop_num_iters_filepath)
        iters = iters[iters["Loop"]==kernel].drop("Loop", axis=1)
        simd_len = 1 if simd_failed else get_attribute(output_dirpath, "SIMD len")

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

            f_nonzero = papi_and_iters["NumIters"] > 0
            papi_and_iters_nonzero = papi_and_iters[papi_and_iters["NumIters"] > 0]
            if papi_and_iters_nonzero.shape[0] == 0:
                ## No rows left
                return -1
            papi_and_iters_nonzero["InsPerIter"] = papi_and_iters_nonzero["Count"] / papi_and_iters_nonzero["NumIters"]
            ins_per_iter = papi_and_iters_nonzero.loc[0,"InsPerIter"]

    if verbose and ins_per_iter == -1:
        print("WARNING: Failed to calculate ins-per-iter for loop '{0}'".format(kernel))

    return ins_per_iter

def get_isa_group(output_dirpath):
    att_filepath = os.path.join(output_dirpath, "Attributes.csv")
    if not os.path.isfile(att_filepath):
        raise Exception("File not found: " + att_filepath)
    attributes = load_attributes(output_dirpath)
    if "ISA group" in attributes:
        return attributes["ISA group"]

    isaGrp = "Intel"
    raw_asm_filepath = os.path.join(output_dirpath, "objects", "flux_loops.o.raw-asm")
    if not os.path.isfile(raw_asm_filepath):
        raw_asm_filepath = os.path.join(output_dirpath, "objects", "flux_vecloops.o.raw-asm")
    for line in open(raw_asm_filepath):
      if "aarch64" in line:
        isaGrp = "ARM"
        break
    return isaGrp

def get_simd_len(output_dirpath):
    run_filepath = os.path.join(output_dirpath, "run.sh")
    if not os.path.isfile(run_filepath):
        raise Exception("Need run.sh file to get attribute '{0}: {1}".format(field, run_filepath))

    value = None
    for line in open(run_filepath):
        res = re.match(".*DBLS_PER_SIMD=([0-9]+).*", line)
        if res:
            value = int(res.groups()[0])
            break
    if value is None:
        raise Exception("Could not determine SIMD length in file {0}".format(run_filepath))

    if value == 0:
        raise Exception("get_simd_len() returning zero")
    return value

def get_simd(output_dirpath):
    run_filepath = os.path.join(output_dirpath, "run.sh")
    if not os.path.isfile(run_filepath):
        raise Exception("Need run.sh file to get attribute '{0}: {1}".format(field, run_filepath))

    value = None
    for line in open(run_filepath):
        res = re.match(".*-DSIMD[ \"].*", line)
        if res:
            value = True
            break
    if value is None:
        raise Exception("Could not determine SIMD status in file {0}".format(run_filepath))
    return value

def get_simd_ca_strategy(output_dirpath):
    run_filepath = os.path.join(output_dirpath, "run.sh")
    if not os.path.isfile(run_filepath):
        raise Exception("Need run.sh file to get attribute '{0}: {1}".format(field, run_filepath))

    value = None
    for line in open(run_filepath):
        res = re.match(".*-DCOLOURED_CONFLICT_AVOIDANCE[ \"].*", line)
        if res:
            value = "colour"
            break

        resm = re.match(".*-DMANUAL[ \"].*", line)
        resmg = re.match(".*-DMANUAL_GATHER[ \"].*", line)
        resms = re.match(".*-DMANUAL_SCATTER[ \"].*", line)
        if resm or (resmg and resms):
            value = "manual"
            break
        elif resmg:
            value = "manual-gather"
            break
        elif resms:
            value = "manual-scatter"
            break
    if value is None:
        raise Exception("Could not determine SIMD status in file {0}".format(run_filepath))
    return value

def get_attribute_from_runsh(output_dirpath, field):
    run_filepath = os.path.join(output_dirpath, "run.sh")
    if not os.path.isfile(run_filepath):
        raise Exception("Need run.sh file to get attribute '{0}: {1}".format(field, run_filepath))

    value = None
    for line in open(run_filepath):
        if "{0}=".format(field) in line:
            value = '='.join(line.split('=')[1:])
            break
    if value is None:
        raise Exception("Could not find attribute '{0}' assignment in file {1}".format(field, run_filepath))
    return value

def get_attribute(output_dirpath, field):
    if field in ["flags_final"]:
        return get_attribute_from_runsh(output_dirpath, field)

    value = None
    att_filepath = os.path.join(output_dirpath, "Attributes.csv")
    if not os.path.isfile(att_filepath):
        ## Try getting values from run.sh script:
        if field == "SIMD len":
            return get_simd_len(output_dirpath)
        elif field == "SIMD":
            return get_simd(output_dirpath)
        elif field == "SIMD conflict avoidance strategy":
            return get_simd_ca_strategy(output_dirpath)

        run_field = ""
        if field == "CC":
            run_field = "compiler"
        if run_field != "":
            return get_attribute_from_runsh(output_dirpath, run_field)

        raise Exception("File not found: " + att_filepath)

    attributes = load_attributes(output_dirpath)
    if field in attributes:
        return attributes[field]

    if value is None:
        raise Exception("Cannot infer field '{0}'".format(field))
    
    return value

def set_loop_attribute(output_dirpath, loop_name, attribute, value):
    att_filepath = os.path.join(output_dirpath, "LoopAttributes.csv")
    if not os.path.isfile(att_filepath):
        row = pd.DataFrame([[loop_name, attribute, value]], columns=["Loop", "Attribute", "Value"])
        row.to_csv(att_filepath, index=False)
        return

    if (attribute == "SIMD failed") and (not isinstance(value, bool)):
        raise Exception("Value for attribute '{0}' should be Boolean, not {1} (value={2})".format(attribute, type(value), value))
    if (attribute == "True SIMD len") and (not isinstance(value, int)):
        raise Exception("Value for attribute '{0}' should be Integer, not {1} (value={2})".format(attribute, type(value), value))

    att_df = clean_pd_read_csv(att_filepath)
    value = str(value)
    f = np.logical_and(att_df["Loop"]==loop_name, att_df["Attribute"]==attribute)
    if sum(f) > 1:
        raise Exception("LoopAttributes.csv contains multiple values for Loop-Attibute pair {0}-{1}".format(loop_name, attribute))
    if sum(f) == 1:
        old_value = att_df.loc[f,"Value"].iloc[0]
        att_df.loc[f,"Value"] = value
    else:
        row = pd.DataFrame([[loop_name, attribute, value]], columns=["Loop", "Attribute", "Value"])
        att_df = att_df.append(row)
    att_df.to_csv(att_filepath, index=False)

def get_loop_attribute(output_dirpath, loop_name, attribute):
    att_filepath = os.path.join(output_dirpath, "LoopAttributes.csv")
    if not os.path.isfile(att_filepath):
        raise Exception("Loop attributes csv does not exist: {0}".format(att_filepath))

    att_df = clean_pd_read_csv(att_filepath)

    f = np.logical_and(att_df["Loop"]==loop_name, att_df["Attribute"]==attribute)
    if sum(f) == 0:
        raise Exception("LoopAttributes.csv contains no value for Loop-Attibute pair {0}-{1}".format(loop_name, attribute))
    if sum(f) > 1:
        raise Exception("LoopAttributes.csv contains multiple values for Loop-Attibute pair {0}-{1}".format(loop_name, attribute))

    value = att_df[f]["Value"].iloc[0]
    value = ast.literal_eval(value)

    return value

def do_logs_report_simd_fail(output_dirpath, func_name, compiler=None):
    if get_attribute(output_dirpath, "SIMD") == "N":
        return False

    valid_funcs = ["compute_flux_edge", "compute_stream", "unstructured_stream", "unstructured_compute"]
    if not func_name in valid_funcs:
        raise Exception("do_logs_report_simd_fail() file prefix '{0}' should be one of: {1}".format(func_name, valid_funcs))
    if func_name == "compute_flux_edge":
        file_prefix = "flux"
    else:
        file_prefix = func_name

    try:
        val = get_loop_attribute(output_dirpath, func_name, "SIMD failed")
        return val
    except:
        pass

    failed = False

    if compiler is None:
        compiler = get_attribute(output_dirpath, "CC")

    if compiler in ["clang", "arm"]:
        log_filepath =  os.path.join(output_dirpath, "objects", "{0}_vecloops.o.log".format(file_prefix))
        if not os.path.isfile(log_filepath):
            log_filepath =  os.path.join(output_dirpath, "objects", "{0}_vecloop.o.log".format(file_prefix))
        f = open(log_filepath, "r")
        print("Parsing log: {0}".format(log_filepath))
        simd_fail_found = False
        simd_success_found = False
        loop_simd_fail = False
        line_num = None
        loop_line_num = None
        for i,line in enumerate(f):
            res = re.match(".*{0}_vecloops\.cpp:([0-9]+):.*".format(file_prefix), line)
            if res:
                line_num = int(res.groups()[0])

                if "loop not vectorized" in line:
                    simd_fail_found = True
                    simd_success_found = False
                    if (not loop_line_num is None) and loop_line_num == line_num:
                        print("'{0}' loop seen previously, and 'not vectorized' on line {1}".format(func_name, i))
                        break
                    ## Source code snippet should be on next line, confirming which loop
                    continue
                elif "vectorized loop" in line:
                    simd_success_found = True
                    simd_fail_found = False
                    if not loop_line_num is None and loop_line_num == line_num:
                        print("'{0}' loop seen previously, and 'vectorized loop' on line {1}".format(func_name, i))
                        break
                    ## Source code snippet should be on next line, confirming which loop
                    continue
                else:
                    ## Discard remarks on vectorization
                    simd_success_found = False
                    simd_fail_found = False

            if "for (long i=flux_loop_start" in line:
                loop_line_num = line_num
                if simd_fail_found:
                    loop_simd_fail = True
                    break
                if simd_success_found:
                    break
            simd_fail_found = False
            simd_success_found = False

        if loop_simd_fail:
            failed = True

    elif compiler == "gnu":
        log_filepath =  os.path.join(output_dirpath, "objects", "gcc-opt-info")
        if os.path.isfile(log_filepath):
            f = open(log_filepath, "r")
            for line in f:
                if re.match("src/Kernels_vectorised/{0}_veckernel.h.*not vectorized".format(file_prefix), line):
                    failed = True
                    break

    elif compiler == "cray":
        log_filepath =  os.path.join(output_dirpath, "objects", "{0}_vecloops.lst".format(file_prefix))
        if not os.path.isfile(log_filepath):
            log_filepath =  os.path.join(output_dirpath, "objects", "{0}_vecloop.lst".format(file_prefix))
        f = open(log_filepath, "r")
        for line in f:
            if "for (long i=flux_loop_start" in line:
                if not re.match(".*V.*for", line):
                    failed = True
                    break

    elif compiler == "fujitsu":
        log_filepath =  os.path.join(output_dirpath, "objects", "{0}_vecloops.lst".format(file_prefix))
        if not os.path.isfile(log_filepath):
            log_filepath =  os.path.join(output_dirpath, "objects", "{0}_vecloop.lst".format(file_prefix))
        f = open(log_filepath, "r")
        ## Sometimes the SIMD indicator appears on line following the for(..):
        for_loop_seen = False
        failed = True ## Assume fail
        for line in f:
            if "for (long i=flux_loop_start" in line:
                for_loop_seen = True
                if re.match(".*v.*for", line):
                    failed = False
                    break
            elif for_loop_seen:
                if re.match(".*v.*{", line):
                    failed = False
                break

    set_loop_attribute(output_dirpath, func_name, "SIMD failed", failed)

    return failed

def analyse_object_files():
    print("Analysing object files")
    global verbose

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

            if verbose:
                print("Processing: " + output_dirpath)

            loop_to_object = {}

            run_filename_candidates = [x+".batch" for x in ["slurm", "pbs", "moab", "lsf"]]
            run_filename_candidates.append("run.sh")
            for rhc in run_filename_candidates:
                if os.path.isfile(os.path.join(output_dirpath, rhc)):
                    run_filename = rhc
                    break
            if run_filename == "":
                raise IOError("Cannot find a run script for run: " + output_dirpath)

            compile_info["compiler"] = get_attribute(output_dirpath, "CC")
            simd_len = get_attribute(output_dirpath, "SIMD len")
            simd = get_attribute(output_dirpath, "SIMD")
            if simd == "N":
                simd_len = 1
                compile_info["SIMD failed"] = False
            compile_info["SIMD len"] = simd_len
            avx512cd_required = compile_info["compiler"] == "intel" and \
                                simd == "Y" and \
                                get_attribute(output_dirpath, "Instruction set") == "AVX512" and \
                                get_attribute(output_dirpath, "SIMD conflict avoidance strategy") == "None"
            compile_info["SIMD CA scheme"] = get_attribute(output_dirpath, "SIMD conflict avoidance strategy")
            compile_info["scatter loop present"] = "manual" in compile_info["SIMD CA scheme"].lower() and "scatter" in compile_info["SIMD CA scheme"].lower()
            compile_info["gather loop present"]  = "manual" in compile_info["SIMD CA scheme"].lower() and "gather"  in compile_info["SIMD CA scheme"].lower()
            flags = get_attribute(output_dirpath, "flags_final")
            NVAR = 5
            NDIM = 3
            compile_info["gather loop numLoads" ] = NVAR*2 + 2 ## 10x variables, 2x node ids
            compile_info["gather loop numStores"] = NVAR*2 ## 10x variables
            compile_info["gather loop numLoads" ] += NDIM  ## 3D edge vector
            compile_info["gather loop numStores"] += NDIM  ## 3D edge vector
            if "FLUX_PRECOMPUTE_EDGE_WEIGHTS" in flags:
                compile_info["gather loop numLoads" ] += 1 ## edge weight
                compile_info["gather loop numStores"] += 1 ## edge weight
            compile_info["scatter loop numLoads" ] = NVAR*2 + NVAR*2 + 2 ## 10x variables, 10x fluxes, 2x node ids
            compile_info["scatter loop numStores"] = NVAR*2              ## 10x fluxes

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

                if simd == "Y":
                    compile_info["SIMD failed"] = do_logs_report_simd_fail(output_dirpath, k_pretty, compile_info["compiler"])
                    compile_info["SIMD len"] = simd_len

                ins_per_iter = calc_ins_per_iter(output_dirpath, k)
                if ins_per_iter == -1:
                    ## Kernel 'k' was not executed
                    continue
                if isnan(ins_per_iter):
                    raise Exception("calc_ins_per_iter(k='{0}') returned NaN, investigate".format(k))

                obj_filepath = os.path.join(output_dirpath, "objects", loop_to_object[k])
                try:
                    loop_asm = extract_loop_kernel_from_obj(obj_filepath, compile_info, ins_per_iter, k, avx512cd_required, 10)
                except:
                    raise
                    # continue
                loop_tally = loop_asm.count_loop_instructions()
                df_data = [ [k_pretty, k, float(v)] for k,v in loop_tally.items()]
                loop_tally_df = pd.DataFrame(df_data, columns=["Loop", "Instruction", "Count"])

                arch = loop_asm.metadata["ARCH"]

                if "SIMD failed" in loop_asm.metadata.keys():
                    simd_failed = loop_asm.metadata["SIMD failed"]
                else:
                    simd_failed = compile_info["SIMD failed"]
                set_loop_attribute(output_dirpath, k_pretty, "SIMD failed", simd_failed)

                if "True SIMD len" in loop_asm.metadata.keys():
                    true_simd_len = loop_asm.metadata["True SIMD len"]
                elif simd_failed:
                    true_simd_len = 1
                else:
                    true_simd_len = compile_info["SIMD len"]
                set_loop_attribute(output_dirpath, k_pretty, "True SIMD len", true_simd_len)

                # if "Unroll factor" in loop_asm.metadata.keys():
                #     unroll_factor = loop_asm.metadata["Unroll factor"]
                # else:
                #     unroll_factor = -1
                # set_loop_attribute(output_dirpath, k_pretty, "Unroll factor", unroll_factor)
                set_loop_attribute(output_dirpath, k_pretty, "#edges per asm pass", loop_asm.metadata["#edges per asm pass"])

                if loops_tally_df is None:
                    loops_tally_df = loop_tally_df
                else:
                    loops_tally_df = loops_tally_df.append(loop_tally_df)

            if not loops_tally_df is None:
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

                target_is_aarch64 = (arch == Architectures.AARCH64)

                loops_tally_categorised_df = None
                loops = Set(loops_tally_df["Loop"])
                for l in loops:
                    df = loops_tally_df[loops_tally_df["Loop"]==l]
                    tally = {}
                    for i in range(df.shape[0]):
                        row = df.iloc[i]
                        tally[row["Instruction"]] = row["Count"]
                    df_metadata_cols = [c for c in df.columns.values if not c in ["Instruction", "Count"]]
                    df_metadata = df[df_metadata_cols].drop_duplicates()

                    if target_is_aarch64:
                        tally_categorised = categorise_aggregated_instructions_tally_dict(tally, is_aarch64=True)
                    else:
                        tally_categorised = categorise_aggregated_instructions_tally_dict(tally, is_intel64=True)

                    ## Don't care about instructions that go straight to load/store execution:
                    for k in ["eu.address", "eu.simd_address"]:
                        if k in tally_categorised:
                            del tally_categorised[k]

                    num_instructions = 0
                    for k in tally_categorised.keys():
                        if k.startswith("eu."):
                            num_instructions += tally_categorised[k]
                    simd_failed = bool(get_loop_attribute(output_dirpath, l, "SIMD failed"))
                    true_simd_len = int(get_loop_attribute(output_dirpath, l, "True SIMD len"))
                    width = int(get_loop_attribute(output_dirpath, l, "#edges per asm pass"))
                    instructions_per_edge = round(num_instructions/width, 1)
                    tally_categorised["#edges per asm pass"] = width
                    tally_categorised["#instructions per edge"] = instructions_per_edge

                    df_data = [ [l, k, float(v)] for k,v in tally_categorised.items()]
                    df = pd.DataFrame(df_data, columns=["Loop", "Instruction", "Count"])
                    df = df.merge(df_metadata)
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
                            if not "ISA group" in df["Attribute"]:
                                isaGrp = get_isa_group(root)
                                df = df.append({"Attribute":"ISA group", "Value":isaGrp}, ignore_index=True)

                        df['Run ID'] = run_id
                        # df = pd.pivot(df, index="Run ID", columns="Attribute", values="Value")

                        new_col_ordering_template = ["Run ID", "Attribute", "Event", "MG level", "Loop", "ThreadNum", "Instruction", "Value", "Count", "Time", "NumIters"]
                        new_col_ordering = [c for c in new_col_ordering_template if c in df.columns.values]
                        if len(new_col_ordering) != len(df.columns.values):
                            raise Exception("cat={0}: New column ordering is missing {1} columns: {2}".format(cat, len(df.columns.values)-len(new_col_ordering), Set(df.columns.values).difference(new_col_ordering)))
                        df = df[new_col_ordering]

                        attributes_cols = ["SIMD failed", "True SIMD len"]
                        attributes_cols_defaults = {"SIMD failed":False, "True SIMD len":-1}
                        if not cat == "Attributes":
                            ## Also merge in loop-attributes:
                            loops_att_df_filepath = os.path.join(root, "LoopAttributes.csv")
                            if os.path.isfile(loops_att_df_filepath):
                                loops_att_df = clean_pd_read_csv(loops_att_df_filepath)
                                loops_att_df = loops_att_df[loops_att_df["Attribute"].isin(attributes_cols)]
                                pivot_id_col = "Attribute" ; pivot_val_col = "Value"
                                index_cols = [c for c in loops_att_df.columns.values if not c in [pivot_id_col, pivot_val_col]]
                                loops_att_df = loops_att_df.pivot(index="Loop", columns="Attribute", values="Value")
                                loops_att_df = loops_att_df.reset_index(drop=False)
                                df = df.merge(loops_att_df, "left")
                                for c in attributes_cols:
                                    df.loc[df[c].isna(),c] = attributes_cols_defaults[c]

                        if agg_dfs[cat] is None:
                            agg_dfs[cat] = df
                        else:
                            for c in attributes_cols:
                                cv = attributes_cols_defaults[c]
                                if (c in agg_dfs[cat].columns.values) and (not c in df.columns.values):
                                    df[c] = cv
                                elif (not c in agg_dfs[cat].columns.values) and (c in df.columns.values):
                                    agg_dfs[cat][c] = cv
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
    if att_pruned_df.shape[0] == 0:
        agg_dfs["Attributes"] = None
    else:
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

        df_grps = df_thread_means.groupby(job_id_colnames)
        df_run_means = df_grps.mean().reset_index()

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

    ## Finally, create wide versions of select key-value tables:
    cat = "instruction-counts-categorised"
    df_in_filepath = os.path.join(prepared_output_dirpath,cat+".csv")
    if os.path.isfile(df_in_filepath):
        pivot_id_col = "Instruction"
        pivot_val_col = "Count"
        print("Reshaping " + cat)
        df_out_filepath = os.path.join(prepared_output_dirpath,cat+".wide.csv")
        df = clean_pd_read_csv(df_in_filepath)
        index_cols = [c for c in df.columns.values if not c in [pivot_id_col, pivot_val_col]]
        df2 = df.pivot_table(index=index_cols, columns=pivot_id_col, values=pivot_val_col)
        df2 = df2.reset_index(drop=False)

        times_filepath = os.path.join(prepared_output_dirpath, "Times.mean.csv")
        if os.path.isfile(times_filepath):
            times_df = pd.read_csv(times_filepath, keep_default_na=False)
            times_df = times_df[times_df["Loop"].isin(df2["Loop"].values)]
            times_df = times_df.drop("MG level", axis=1)
            job_id_colnames = get_job_id_colnames(times_df)
            data_colnames = [c for c in times_df.columns.values if not c in job_id_colnames]
            times_df = times_df.groupby(job_id_colnames).sum().reset_index()
            df2 = df2.merge(times_df, validate="one_to_one")
            # Reorder columns:
            final_cols = ["#instructions per edge", "Time"]
            reserved_cols = ["#edges per asm pass"]
            insn_cols = [c for c in df2.columns.values if "eu." in c or "mem." in c]
            id_cols = [c for c in df2.columns.values if not (c in insn_cols+final_cols+reserved_cols)]
            if "#edges per asm pass" in df2.columns.values:
                id_cols.append("#edges per asm pass")
            df2 = df2[id_cols + insn_cols + final_cols]
        df2.to_csv(df_out_filepath, index=False)

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
        if "True SIMD len" in fp_df_sum.columns.values:
            f = fp_df_sum["True SIMD len"] != -1
            fp_df_sum.loc[f,"SIMD len"] = fp_df_sum.loc[f,"True SIMD len"]
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
        iters_df = iters_df[iters_df["NumIters"] != 0.0]
        iters_df = safe_pd_filter(iters_df, "Metric", "#iterations_SUM")
        if "SIMD len" in iters_df.columns.values:
            ## Need the number of SIMD iterations. Currently 'iterations_SUM'
            ## counts the number of serial iterations.
            if "True SIMD len" in iters_df.columns.values:
                f = iters_df["True SIMD len"] != -1
                iters_df.loc[f,"SIMD len"] = iters_df.loc[f,"True SIMD len"]
            simd_mask = np.invert(iters_df["SIMD failed"])
            data_colnames = get_data_colnames(iters_df)
            for dc in data_colnames:
                iters_df.loc[simd_mask, dc] /= iters_df.loc[simd_mask, "SIMD len"]
        flux_iters_df = iters_df[iters_df["Loop"] == "compute_flux_edge"]

        fp_counts = count_fp_ins()
        if not fp_counts is None:
            ## Calculate GFLOPs
            flops_total = flux_iters_df.merge(fp_counts.drop("FP ins/iter", axis=1))
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
                fp_total = flux_iters_df.merge(fp_counts.drop("FLOPs/iter", axis=1))
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
        gbsec_data = gbsec_data[gbsec_data["Value"] != 0.0]
        data_all = safe_pd_append(data_all, gbsec_data)
    if "Runtime" in metrics and "GFLOPs" in metrics:
        ## Calculate GFLOPs/sec
        runtime_data = data_all[data_all["Metric"]=="Runtime"]
        gflops_data = data_all[data_all["Metric"]=="GFLOPs"]
        gflops_data = safe_frame_divide(gflops_data, runtime_data.drop("Metric", axis=1))
        gflops_data["Metric"] = "GFLOPs/sec"
        gflops_data = gflops_data[gflops_data["Value"] != 0.0]
        data_all = safe_pd_append(data_all, gflops_data)
    if "GFLOPs" in metrics and "dram read GB" in metrics:
        ## Calculate Flops/Byte
        gflops_data = data_all[data_all["Metric"]=="GFLOPs"]
        gb_data = data_all[data_all["Metric"]=="dram read GB"]
        arith_intensity_data = safe_frame_divide(gflops_data, gb_data.drop("Metric", axis=1))
        arith_intensity_data["Metric"] = "Flops/Byte"
        arith_intensity_data = arith_intensity_data[arith_intensity_data["Value"] != 0.0]
        data_all = safe_pd_append(data_all, arith_intensity_data)
    if "Runtime" in metrics and "PAPI_TOT_CYC.THREADS_MAX" in metrics:
        ## Calculate GHz
        cyc_data = data_all[data_all["Metric"]=="PAPI_TOT_CYC.THREADS_MAX"]
        runtime_data = data_all[data_all["Metric"]=="Runtime"]
        ghz = safe_frame_divide(cyc_data, runtime_data.drop("Metric", axis=1))
        ghz["Value"] /= 1e9
        ghz["Metric"] = "GHz"
        ghz = ghz[ghz["Value"] != 0.0]
        data_all = safe_pd_append(data_all, ghz)
    if "GFLOPs" in metrics and "PAPI_TOT_CYC.THREADS_MAX" in metrics:
        ## Calculate flops/cycle
        # cyc_data = data_all[data_all["metric"]=="PAPI_TOT_CYC.THREADS_MAX"]
        cyc_data = data_all[data_all["Metric"]=="PAPI_TOT_CYC"]
        gflops_data = data_all[data_all["Metric"]=="GFLOPs"]
        flops_per_cyc = safe_frame_divide(gflops_data, cyc_data.drop("Metric", axis=1))
        flops_per_cyc["Value"] *= 1e9
        flops_per_cyc["Metric"] = "Flops/Cycle"
        flops_per_cyc = flops_per_cyc[flops_per_cyc["Value"] != 0.0]
        data_all = safe_pd_append(data_all, flops_per_cyc)

    # # Drop runtime:
    # f = data_all["Metric"] == "Runtime"
    # data_all = data_all[np.logical_not(f)]

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
