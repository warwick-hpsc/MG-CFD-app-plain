import os, shutil, stat, sys
import json, argparse, re
import itertools
import math
from tempfile import NamedTemporaryFile
import imp

import sys
pyv = sys.version_info[0]

script_dirpath = os.path.join(os.getcwd(), os.path.dirname(__file__))
template_dirpath = os.path.join(os.path.dirname(script_dirpath), "run-templates")
app_dirpath = os.path.dirname(script_dirpath)

imp.load_source('utils', os.path.join(script_dirpath, "utils.py"))
from utils import *

js_to_filename = {}
js_to_filename[""] = ""
js_to_filename["slurm"] = "slurm.sh"
js_to_filename["moab"] = "moab.sh"
js_to_filename["lsf"] = "lsf.sh"
js_to_filename["pbs"] = "pbs.sh"

js_to_submit_cmd = {}
js_to_submit_cmd[""] = ""
js_to_submit_cmd["slurm"] = "sbatch"
js_to_submit_cmd["moab"] = "msub"
js_to_submit_cmd["lsf"] = "bsub"
js_to_submit_cmd["pbs"] = "qsub -V"

defaults = {}
# Compilation:
defaults["cpp wrapper"] = ""
defaults["debug"] = False
defaults["insn set"] = "Host"
defaults["base flags"] = "-DTIME"
defaults["flux flags"] = [""]
defaults["precise fp"] = True
defaults["compile only"] = False
# Job scheduling:
defaults["unit walltime"] = 0.0
defaults["budget code"] = ""
defaults["partition"] = ""
defaults["single batch"] = False
# MG-CFD execution:
defaults["num threads"] = 1
defaults["num repeats"] = 1
defaults["mg cycles"] = 10
defaults["validate result"] = False
defaults["measure mem bound"] = False
defaults["measure compute bound"] = False
defaults["run synthetic compute"] = False
defaults["min mesh multi"] = 1
# Optimisation:
defaults["simd mode"] = False
defaults["simd CA scheme"] = ""
defaults["renumber"] = False

def get_key_value(profile, cat, key, ensure_list=False):
    v = None
    if cat in profile and key in profile[cat]:
        v = profile[cat][key]
    else:
        if key in defaults:
            v = defaults[key]
        else:
            raise Exception("Mandatory key '{0}' not present in cat '{1}' of json".format(key, cat))

    if ensure_list and (not v is None) and (not isinstance(v, list)):
        v = [v]

    return v

def py_sed(filepath, from_rgx, to_rgx, delete_if_null=False):
    with open(filepath, "r") as f:
        lines = f.readlines()
    with open(filepath, "w") as f:
        for line in lines:
            if re.search(from_rgx, line):
                if to_rgx is None or to_rgx == "" and delete_if_null:
                    ## Discard this line
                    pass
                else:
                    if isinstance(to_rgx, str):
                        f.write(re.sub(from_rgx,     to_rgx,  line))
                    else:
                        f.write(re.sub(from_rgx, str(to_rgx), line))
            else:
                f.write(line)

def py_grep(filepath, str):
    with open(filepath, "r") as f:
        lines = f.readlines()
    for line in lines:
        if line.find(str) != -1:
            return True

    return False

def delete_folder_contents(dirpath):
    print("Deleting contents of folder: " + dirpath)
    for f in os.listdir(dirpath):
        fp = os.path.join(dirpath, f)
        if os.path.isdir(fp):
            shutil.rmtree(fp)
        else:
            os.remove(fp)

def prune_flux_flags_permutations(flux_flags_permutations):
    flux_flags_permutations_pruned = []
    for perm in flux_flags_permutations:
        other_flux_options_present = False
        for i in perm:
            if "FLUX_" in i:
                other_flux_options_present = True
                break
        perm_clean = []
        for i in perm:
            if i == "" and (i in perm_clean or other_flux_options_present):
                continue
            else:
                perm_clean.append(i)
        if len(perm_clean) > 0:
            flux_flags_permutations_pruned.append(perm_clean)

    # Remove duplicates:
    flux_flags_permutations_pruned.sort()
    flux_flags_permutations_pruned = list(k for k,_ in itertools.groupby(flux_flags_permutations_pruned))

    if len(flux_flags_permutations_pruned) == 0:
        flux_flags_permutations_pruned = [""]

    return flux_flags_permutations_pruned

def powerset(iterable):
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))

def concat_compile_flags(flags):
    if isinstance(flags, str) or (pyv==2 and isinstance(flags, unicode)):
        flags = flags.split(' ')

    s = ""
    for f in flags:
        if f != "":
            if not f.startswith("-D"):
                f2 = "-D"+f
            else:
                f2 = f
            if s == "":
                s = f2
            else:
                s += " " + f2
    return s

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--json', required=True)
    args = parser.parse_args()
    with open(args.json, 'r') as f:
        profile = json.load(f)

    ## Read file/folder paths and prepare folders:
    jobs_dir = profile["setup"]["jobs dir"]
    if jobs_dir[0] != '/':
        jobs_dir = os.path.join(os.getcwd(), jobs_dir)
    if not os.path.isdir(jobs_dir):
        os.mkdir(jobs_dir)
    else:
        delete_folder_contents(jobs_dir)

    data_dirpath = profile["run"]["data dirpath"]
    if data_dirpath[0] != '/':
        data_dirpath = os.path.join(app_dirpath, data_dirpath)

    ## Read parameters from json:
    job_queue = get_key_value(profile, "setup", "partition")
    budget_code = get_key_value(profile, "setup", "budget code")
    js = get_key_value(profile, "setup", "job scheduler")
    single_batch = get_key_value(profile, "setup", "single batch") 

    compilers = get_key_value(profile, "compile", "compiler", True)
    cpp_wrapper = get_key_value(profile, "compile", "cpp wrapper")
    compile_only = get_key_value(profile, "compile", "compile only")
    debug = get_key_value(profile, "compile", "debug")
    if "insn sets" in profile["compile"].keys():
        ## Backwards compatibility
        isas = get_key_value(profile, "compile", "insn sets", True)
    else:
        isas = get_key_value(profile, "compile", "insn set", True)
    base_flags = get_key_value(profile, "compile", "base flags", True)
    precise_fp = get_key_value(profile, "optimisation", "precise fp", True)

    if "permutable flags" in profile["compile"]:
        ## Backwards compatibility
        perm_flags = get_key_value(profile, "compile", "permutable flags")
        flux_flags_permutations = prune_flux_flags_permutations(itertools.product(*perm_flags))
    else:
        flux_flags = get_key_value(profile, "compile", "flux flags", True)
        flux_flags_permutations = prune_flux_flags_permutations(powerset(flux_flags))

    threads = get_key_value(profile, "run", "num threads", True)
    num_repeats = get_key_value(profile, "run", "num repeats")
    mg_cycles = get_key_value(profile, "run", "mg cycles")
    validate = get_key_value(profile, "run", "validate result")
    renumber_modes = get_key_value(profile, "optimisation", "renumber", True)
    measure_mem_bound = get_key_value(profile, "run", "measure mem bound")
    measure_compute_bound = get_key_value(profile, "run", "measure compute bound")
    run_synthetic_compute = get_key_value(profile, "run", "run synthetic compute")
    mgcfd_unit_runtime_secs = get_key_value(profile, "run", "unit walltime")
    min_mesh_multi = get_key_value(profile, "run", "min mesh multi")

    simd_modes = get_key_value(profile, "optimisation", "simd mode", True)
    if True in simd_modes:
        simd_ca_schemes = get_key_value(profile, "optimisation", "simd CA scheme", True)
        simd_lens = get_key_value(profile, "optimisation", "simd len", True)
    else:
        simd_ca_schemes = None
        simd_lens = None

    iteration_space = {}
    if not compilers is None:
        iteration_space["compiler"] = compilers
    if not isas is None:
        iteration_space["isa"] = isas
    if not base_flags is None:
        iteration_space["base flags"] = base_flags
    if not precise_fp is None:
        iteration_space["precise fp"] = precise_fp
    if not flux_flags_permutations is None:
        iteration_space["flux flags"] = flux_flags_permutations
    if not threads is None:
        iteration_space["threads"] = threads
    if not renumber_modes is None:
        iteration_space["renumber"] = renumber_modes
    if not simd_modes is None:
        iteration_space["simd mode"] = simd_modes
    if not simd_ca_schemes is None:
        iteration_space["simd CA scheme"] = simd_ca_schemes
    if not simd_lens is None:
        iteration_space["simd len"] = simd_lens

    # with open(os.path.join(jobs_dir, "papi.conf"), "w") as f:
    #     if "papi events" in profile["run"].keys():
    #         for e in profile["run"]["papi events"]:
    #             f.write("{0}\n".format(e))
    #     else:
    #         # Default events:
    #         f.write("PAPI_TOT_INS\n")
    #         f.write("PAPI_TOT_CYC\n")
    ## Use 'papi_event_chooser' to check for incompatibility between events
    papi_preset_events = []
    papi_native_events = []
    papi_requested_events = []
    batched_papi_events = []
    if "papi events" in profile["run"].keys():
        events = profile["run"]["papi events"]
        if events[0] is list:
            ## User has already grouped PAPI events into batches:
            batched_papi_events = events
        for e in events:
            if e.startswith("PAPI_"):
                papi_preset_events.append(e)
            else:
                papi_native_events.append(e)
    else:
        papi_preset_events = ["PAPI_TOT_INS", "PAPI_TOT_CYC"]
    batched_papi_events = batch_papi_events(papi_preset_events, papi_native_events)
    if len(batched_papi_events) > 0:
        iteration_space["batched papi events"] = batched_papi_events

    iterables = itertools.product(*iteration_space.values())
    iterables_labelled = []
    iteration_keys = list(iteration_space.keys())
    for item in iterables:
        item_dict = {}
        for i in range(len(item)):
          item_dict[iteration_keys[i]] = item[i]
        iterables_labelled.append(item_dict)
    ## Process incompatible options:
    iterables_labelled_processed = []
    for item in iterables_labelled:
        if not item["simd mode"]:
            item["simd len"] = 1
            item["simd CA scheme"] = ""
        else:
            isa = item["isa"]
            simd_len = item["simd len"]
            if isa.startswith("SSE42"):
                if simd_len > 2:
                    item["simd len"] = 2
            elif isa.startswith("AVX") and not isa == "AVX512":
                if simd_len > 4:
                    item["simd len"] = 4

        if item["compiler"] == "intel" and item["precise fp"]:
            disable_precise_fp = False
            if item["isa"] == "AVX512" or item["isa"] == "KNL-AVX512":
                disable_precise_fp = True
            elif item["isa"] == "Host" and ("avx512cd" in subprocess.check_output(['lscpu'])):
                disable_precise_fp = True
            if disable_precise_fp:
                print("WARNING: Intel compiler segfaults when compiling precise FP with target AVX512, so disabling precise FP")
                item["precise fp"] = False
        iterables_labelled_processed.append(item)
    ## Prune duplicates:
    iterables_labelled_pruned = []
    for item in iterables_labelled_processed:
        if not item in iterables_labelled_pruned:
            iterables_labelled_pruned.append(item)
    iterables_labelled = iterables_labelled_pruned

    num_jobs = len(iterables_labelled) * num_repeats

    # with open(os.path.join(jobs_dir, "papi.conf"), "w") as f:
    #     if "papi events" in profile["run"].keys():
    #         for e in profile["run"]["papi events"]:
    #             f.write("{0}\n".format(e))
    #     else:
    #         # Default events:
    #         f.write("PAPI_TOT_INS\n")
    #         f.write("PAPI_TOT_CYC\n")
    # Use 'papi_event_chooser' to check for incompatibility between events
    papi_preset_events = []
    papi_native_events = []
    papi_requested_events = []
    batched_papi_events = []
    if "papi events" in profile["run"].keys():
        events = profile["run"]["papi events"]
        if events[0] is list:
            ## User has already grouped PAPI events into batches:
            batched_papi_events = events
        for e in events:
            # papi_events.append(e)
            if e.startswith("PAPI_"):
                papi_preset_events.append(e)
            else:
                papi_native_events.append(e)
    else:
        papi_preset_events = ["PAPI_TOT_INS", "PAPI_TOT_CYC"]
    batched_papi_events = batch_papi_events(papi_preset_events, papi_native_events)
    if len(batched_papi_events) > 0:
        iteration_space["batched papi events"] = batched_papi_events

    if js != "":
        js_filename = js_to_filename[js]
    if "batch header filepath" in profile["setup"].keys():
        src_js_filepath = profile["setup"]["batch header filepath"]
        if not os.path.isfile(src_js_filepath):
            raise Exception("Cannot find provided 'batch header filepath': '{0}'".format(src_js_filepath))
    else:
        src_js_filepath = None
    if (src_js_filepath is None) and (js != ""):
        ## Use bundled scheduler options
        src_js_filepath = os.path.join(template_dirpath, js_filename)

    if js == "":
        single_batch = False

    if (not src_js_filepath is None) and py_grep(src_js_filepath, "<COMPILER>") and len(compilers) > 1:
        raise Exception("Single batch is incompatible with multiple compilers")

    ## Init the master job submission script:
    submit_all_filepath = os.path.join(jobs_dir, "submit_all.sh")
    submit_all_file = open(submit_all_filepath, "w")
    submit_all_file.write("#!/bin/bash\n")
    if single_batch:
        run_all_filepath = os.path.join(jobs_dir, "run_all.sh")
        run_all_file = open(run_all_filepath, "w")
        with open(src_js_filepath, "r") as f_in:
            for line in f_in.readlines():
                run_all_file.write(line)
        run_all_file.write("submit_cmd=\"\"")
        run_all_file.write("\n")
        run_all_file.write("set -e\n")
        run_all_file.write("set -u\n")
        run_all_file.write("\n")

        submit_all_file.write("\n")
        submit_all_file.write("{0} {1}".format(js_to_submit_cmd[js], run_all_filepath))

        estimated_total_runtime_secs = 0.0
        nt_max = 1
        ## No longer need to know which job scheduler:
        js = ""
        run_all_file.write("num_jobs={0}\n\n".format(num_jobs))
    else:
        submit_all_file.write("set -e\n")
        submit_all_file.write("set -u\n")
        submit_all_file.write("\n")
        submit_all_file.write("# {0}:\n".format(js))
        submit_all_file.write("submit_cmd=\"{0}\"\n\n".format(js_to_submit_cmd[js]))
        submit_all_file.write("num_jobs={0}\n\n".format(num_jobs))

    n = 0
    for repeat in range(num_repeats):
        for item in iterables_labelled:
            n += 1
            job_id = str(n).zfill(3)
            print("Creating job {0}/{1}".format(n, num_jobs))

            job_dir = os.path.join(jobs_dir, job_id)
            if not os.path.isdir(job_dir):
                os.mkdir(job_dir)

            ## Prepare compilation flags:
            build_flags = concat_compile_flags(item.get("base flags", None))
            p = item.get("flux flags", None)
            flux_flags = concat_compile_flags(p)
            if flux_flags != "":
                build_flags += " " + flux_flags

            renumber = item.get("renumber")

            compiler = item.get("compiler")
            isa = item.get("isa")
            simd = item.get("simd mode")
            if simd:
                build_flags += " -DSIMD"
                build_flags += " -DDBLS_PER_SIMD={0}".format(item.get("simd len"))
                ca_scheme = item.get("simd CA scheme")
                if ca_scheme == "manual":
                    build_flags += " -DMANUAL_SCATTER -DMANUAL_GATHER"
                elif ca_scheme == "manual scatter":
                    build_flags += " -DMANUAL_SCATTER"
                elif ca_scheme == "colour":
                    build_flags += " -DCOLOURED_CONFLICT_AVOIDANCE"
                    build_flags += " -DBIN_COLOURED_VECTORS"

            precise_fp = item.get("precise fp")
            if precise_fp:
                build_flags += " -DPRECISE_FP"

            bin_filename = "euler3d_cpu_double_" + compiler
            bin_filename += build_flags.replace(' ', '')+"-DINSN_SET="+isa+".b"
            bin_filepath = os.path.join(app_dirpath, "bin", bin_filename)

            nt = item.get("threads")

            ## Link to papi config file:
            papi_dest_filepath = os.path.join(job_dir, "papi.conf")
            if ("batched papi events" in item.keys()) and len(item["batched papi events"]) > 0:
                with open(papi_dest_filepath, "w") as f_out:
                    for e in item["batched papi events"]:
                        f_out.write(e+"\n")

            ## Instantiate MG-CFD run script:
            job_run_filepath = os.path.join(job_dir, "run-mgcfd.sh")
            shutil.copyfile(os.path.join(template_dirpath, "run-mgcfd.sh"), job_run_filepath)

            ## Instantiate job scheduling header:
            if not src_js_filepath is None:
                out_js_filepath = os.path.join(job_dir, js_filename)
                shutil.copyfile(src_js_filepath, out_js_filepath)

            ## Combine into a batch run script:
            if js == "":
                batch_filename = "run.sh"
            else:
                batch_filename = js+".batch"
            batch_filepath = os.path.join(job_dir, batch_filename)
            with open(batch_filepath, "w") as f_out:
                if js != "":
                    with open(out_js_filepath, "r") as f_in:
                        for line in f_in.readlines():
                            f_out.write(line)
                    os.remove(out_js_filepath)
                else:
                    f_out.write("#!/bin/bash\n")
                f_out.write("\n\n")
                with open(job_run_filepath, "r") as f_in:
                    for line in f_in.readlines():
                        f_out.write(line)
                os.remove(job_run_filepath)

            ## Now replace variables in script:

            ## - File/dir paths:
            py_sed(batch_filepath, "<RUN_OUTDIR>", job_dir)
            py_sed(batch_filepath, "<APP_DIRPATH>", app_dirpath)
            py_sed(batch_filepath, "<DATA_DIRPATH>", data_dirpath)

            ## - Scheduling:
            py_sed(batch_filepath, "<RUN ID>", job_id)
            py_sed(batch_filepath, "<PARTITION>", job_queue)
            py_sed(batch_filepath, "<RUN_DIR>", job_dir)
            py_sed(batch_filepath, "<BUDGET CODE>", budget_code, True)

            ## - Parallelism:
            py_sed(batch_filepath, "<NUM_THREADS>", nt)
            mesh_multi = min_mesh_multi
            if not "-DFLUX_FISSION" in build_flags:
                ## Duplicate mesh to ensure each thread has a whole copy:
                while (mesh_multi % nt) > 0:
                    mesh_multi += 1
            py_sed(batch_filepath, "<MESH_MULTI>", mesh_multi)
            if single_batch:
                nt_max = max(nt_max, nt)

            ## - Compilation:
            py_sed(batch_filepath, "<COMPILER>", compiler)
            py_sed(batch_filepath, "<ISA>", isa)
            py_sed(batch_filepath, "<BUILD_FLAGS>", build_flags)
            py_sed(batch_filepath, "<CPP_WRAPPER>", cpp_wrapper)
            py_sed(batch_filepath, "<DEBUG>", str(debug).lower())

            ## - Execution:
            py_sed(batch_filepath, "<MG_CYCLES>", mg_cycles)
            py_sed(batch_filepath, "<VALIDATE_RESULT>", str(validate).lower())
            py_sed(batch_filepath, "<RENUMBER>", str(renumber).lower())
            py_sed(batch_filepath, "<MEASURE_MEM_BOUND>", str(measure_mem_bound).lower())
            py_sed(batch_filepath, "<MEASURE_COMPUTE_BOUND>", str(measure_compute_bound).lower())
            py_sed(batch_filepath, "<RUN_SYNTHETIC_COMPUTE>", str(run_synthetic_compute).lower())

            ## - Walltime estimation:
            if single_batch:
                if mgcfd_unit_runtime_secs == 0.0:
                    estimated_total_runtime_secs += 60.0 * 10
                else:
                    estimated_total_runtime_secs += float(mgcfd_unit_runtime_secs*mg_cycles*mesh_multi) / math.sqrt(float(nt))
            else:
                if mgcfd_unit_runtime_secs == 0.0:
                    est_runtime_hours = 0
                    est_runtime_minutes = 10
                else:
                    est_runtime_secs = float(mgcfd_unit_runtime_secs*mg_cycles*mesh_multi) / math.sqrt(float(nt))
                    est_runtime_secs *= 1.2 ## Allow for estimation error
                    est_runtime_secs += 20  ## Add time for file load
                    est_runtime_secs += 60
                    est_runtime_secs = int(round(est_runtime_secs))
                    if pyv == 3:
                        est_runtime_hours = est_runtime_secs//60//60
                        est_runtime_secs -= est_runtime_hours*60*60
                        est_runtime_minutes = est_runtime_secs//60
                        est_runtime_secs -= est_runtime_minutes*60
                    else:
                        est_runtime_hours = est_runtime_secs/60/60
                        est_runtime_secs -= est_runtime_hours*60*60
                        est_runtime_minutes = est_runtime_secs/60
                        est_runtime_secs -= est_runtime_minutes*60
                    if est_runtime_secs > 0:
                        est_runtime_minutes += 1
                        est_runtime_secs = 0
                py_sed(batch_filepath, "<HOURS>", str(est_runtime_hours).zfill(2))
                py_sed(batch_filepath, "<MINUTES>", str(est_runtime_minutes).zfill(2))

            ## Make batch script executable (chmod 755):
            os.chmod(batch_filepath, stat.S_IRWXU | stat.S_IRGRP|stat.S_IXGRP | stat.S_IROTH|stat.S_IXOTH)

            ## Add an entry to submit_all.sh:
            # Copy template 'submit.sh' to a temp file:
            submit_tmp = NamedTemporaryFile(prefix='myprefix')
            if not os.path.isfile(submit_tmp.name):
                print("NamedTemporaryFile() failed to actually create the file")
                sys.exit(-1)
            submit_tmp_filepath = submit_tmp.name
            with open(os.path.join(template_dirpath, "submit.sh"), 'r+b') as f:
                shutil.copyfileobj(f, submit_tmp)
            submit_tmp.seek(0)
            # Instantiate and append to submit_all.sh: 
            py_sed(submit_tmp_filepath, "<RUN_DIRPATH>",    job_dir)
            py_sed(submit_tmp_filepath, "<BATCH_FILENAME>", batch_filename)
            py_sed(submit_tmp_filepath, "<BIN_FILEPATH>",   bin_filepath)
            py_sed(submit_tmp_filepath, "<COMPILE_ONLY>", str(compile_only).lower())
            if not single_batch:
                submit_all_file.write("\n\n")
            with open(submit_tmp_filepath, 'r') as f:
                for line in f:
                    if single_batch:
                        run_all_file.write(line)
                    else:
                        submit_all_file.write(line)
            # Now close (and delete) submit_tmp file
            submit_tmp.close()
            if os.path.isfile(submit_tmp_filepath):
                os.unlink(submit_tmp_filepath)

    submit_all_file.write("\n\n")
    submit_all_file.write('echo "ALL JOBS HAVE BEEN SUBMITTED"\n')

    submit_all_file.close()
    st = os.stat(submit_all_filepath)
    os.chmod(submit_all_filepath, st.st_mode | stat.S_IEXEC)
    if single_batch:
        run_all_file.close()
        st = os.stat(run_all_filepath)
        os.chmod(run_all_filepath, st.st_mode | stat.S_IEXEC)

        ## Fill in scheduling fields:
        est_runtime_secs = estimated_total_runtime_secs
        est_runtime_secs *= 1.2 ## Allow for estimation error
        est_runtime_secs += 20  ## Add time for file load
        est_runtime_secs += 60
        est_runtime_secs = int(round(est_runtime_secs))
        if pyv == 3:
            est_runtime_hours = est_runtime_secs/60/60
            est_runtime_secs -= est_runtime_hours*60*60
            est_runtime_minutes = est_runtime_secs/60
            est_runtime_secs -= est_runtime_minutes*60
        else:
            est_runtime_hours = est_runtime_secs/60/60
            est_runtime_secs -= est_runtime_hours*60*60
            est_runtime_minutes = est_runtime_secs/60
            est_runtime_secs -= est_runtime_minutes*60
        if est_runtime_secs > 0:
            est_runtime_minutes += 1
            est_runtime_secs = 0
        py_sed(run_all_filepath, "<HOURS>", str(est_runtime_hours).zfill(2))
        py_sed(run_all_filepath, "<MINUTES>", str(est_runtime_minutes).zfill(2))

        py_sed(run_all_filepath, "<RUN ID>", 1)
        py_sed(run_all_filepath, "<NUM_THREADS>", nt_max)
        py_sed(run_all_filepath, "<COMPILER>", compilers[0])
        py_sed(run_all_filepath, "<PARTITION>", job_queue)
        py_sed(run_all_filepath, "<BUDGET CODE>", budget_code, True)

    ## Create a little script for listing jobs that failed during execution:
    error_scan_filepath = os.path.join(jobs_dir, "list_errored_jobs.sh")
    error_scan_file = open(error_scan_filepath, "w")
    error_scan_file.write("#!/bin/bash\n")
    error_scan_file.write("find . -name pbs.stderr | while read F ; do wc -l \"$F\" ; done | grep -v \"^0 \"\n")
    error_scan_file.close()
    st = os.stat(error_scan_filepath)
    os.chmod(error_scan_filepath, st.st_mode | stat.S_IEXEC)
