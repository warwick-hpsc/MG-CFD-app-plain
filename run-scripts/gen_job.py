import os, shutil, stat, sys
import json, argparse, re
import itertools
import math
from tempfile import NamedTemporaryFile

script_dirpath = os.path.join(os.getcwd(), os.path.dirname(__file__))
template_dirpath = os.path.join(os.path.dirname(script_dirpath), "run-templates")
app_dirpath = os.path.dirname(script_dirpath)

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
# Job scheduling:
defaults["unit walltime"] = 0.0
defaults["budget code"] = ""
defaults["partition"] = ""
# MG-CFD execution:
defaults["num threads"] = 1
defaults["num repeats"] = 1
defaults["mg cycles"] = 10
defaults["validate result"] = False
defaults["min mesh multi"] = 1

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
    for p in flux_flags_permutations:
        flux_cripple_option_present = False
        other_flux_options_present = False
        for i in p:
            if "FLUX_CRIPPLE" in i:
                flux_cripple_option_present = True
            elif "FLUX_" in i:
                other_flux_options_present = True
        if flux_cripple_option_present and other_flux_options_present:
            ## Discard this permutation, as the FLUX_CRIPPLE option should be
            ## applied without other flux options.
            pass
        else:
            p_clean = []
            for i in p:
                if i == "" and ("" in p_clean or other_flux_options_present or flux_cripple_option_present):
                    continue
                else:
                    p_clean.append(i)
            if len(p_clean) > 0:
                flux_flags_permutations_pruned.append(p_clean)

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
    if isinstance(flags, str) or isinstance(flags, unicode):
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

    compilers = get_key_value(profile, "compile", "compiler", True)
    cpp_wrapper = get_key_value(profile, "compile", "cpp wrapper")
    debug = get_key_value(profile, "compile", "debug")
    if "insn sets" in profile["compile"].keys():
        ## Backwards compatibility
        isas = get_key_value(profile, "compile", "insn sets", True)
    else:
        isas = get_key_value(profile, "compile", "insn set", True)
    base_flags = get_key_value(profile, "compile", "base flags", True)

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
    mgcfd_unit_runtime_secs = get_key_value(profile, "run", "unit walltime")
    min_mesh_multi = get_key_value(profile, "run", "min mesh multi")

    iteration_space = {}
    if not compilers is None:
        iteration_space["compiler"] = compilers
    if not isas is None:
        iteration_space["isa"] = isas
    if not base_flags is None:
        iteration_space["base flags"] = base_flags
    if not flux_flags_permutations is None:
        iteration_space["flux flags"] = flux_flags_permutations
    if not threads is None:
        iteration_space["threads"] = threads

    iterables = itertools.product(*iteration_space.values())
    iterables_labelled = []
    for item in iterables:
        item_dict = {}
        for i in xrange(len(item)):
          item_dict[iteration_space.keys()[i]] = item[i]
        iterables_labelled.append(item_dict)
    ## Prune duplicates:
    iterables_labelled_pruned = []
    for item in iterables_labelled:
        if not item in iterables_labelled_pruned:
            iterables_labelled_pruned.append(item)
    iterables_labelled = iterables_labelled_pruned

    num_jobs = len(iterables_labelled) * num_repeats

    ## Init the master job submission script:
    submit_all_filepath = os.path.join(jobs_dir, "submit_all.sh")
    submit_all_file = open(submit_all_filepath, "w")
    submit_all_file.write("#!/bin/bash\n")
    submit_all_file.write("set -e\n")
    submit_all_file.write("set -u\n")
    submit_all_file.write("\n")
    submit_all_file.write("# {0}:\n".format(js))
    submit_all_file.write("submit_cmd=\"{0}\"\n\n".format(js_to_submit_cmd[js]))
    submit_all_file.write("num_jobs={0}\n\n".format(num_jobs))

    if js != "":
        js_filename = js_to_filename[js]
    if "batch header filepath" in profile["setup"].keys():
        src_js_filepath = profile["setup"]["batch header filepath"]
        if not os.path.isfile(src_js_filepath):
            raise Exception("Cannot find provided 'batch header filepath': '{0}'".format(src_js_filepath))
    else:
        src_js_filepath = None

    with open(os.path.join(jobs_dir, "papi.conf"), "w") as f:
        f.write("PAPI_TOT_INS\n")
        f.write("PAPI_TOT_CYC\n")

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

            compiler = item.get("compiler")
            isa = item.get("isa")

            bin_filename = "euler3d_cpu_double_" + compiler
            bin_filename += build_flags.replace(' ', '')+"-DINSN_SET="+isa+".b"
            bin_filepath = os.path.join(app_dirpath, "bin", bin_filename)

            nt = item.get("threads")

            ## Link to papi config file:
            papi_dest_filepath = os.path.join(job_dir, "papi.conf")
            if os.path.isfile(papi_dest_filepath):
                os.remove(papi_dest_filepath)
            os.symlink(os.path.join(jobs_dir, "papi.conf"), papi_dest_filepath)

            ## Instantiate MG-CFD run script:
            job_run_filepath = os.path.join(job_dir, "run-mgcfd.sh")
            shutil.copyfile(os.path.join(template_dirpath, "run-mgcfd.sh"), job_run_filepath)

            ## Instantiate job scheduling header:
            if (src_js_filepath is None) and (js != ""):
                ## Use bundled scheduler options
                src_js_filepath = os.path.join(template_dirpath, js_filename)
            if not src_js_filepath is None:
                out_js_filepath = os.path.join(job_dir, js_filename)
                shutil.copyfile(src_js_filepath, out_js_filepath)

            ## Combine into a batch submission script:
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

            ## - Compilation:
            py_sed(batch_filepath, "<COMPILER>", compiler)
            py_sed(batch_filepath, "<ISA>", isa)
            py_sed(batch_filepath, "<BUILD_FLAGS>", build_flags)
            py_sed(batch_filepath, "<CPP_WRAPPER>", cpp_wrapper)
            py_sed(batch_filepath, "<DEBUG>", str(debug).lower())

            ## - Execution:
            py_sed(batch_filepath, "<MG_CYCLES>", mg_cycles)
            py_sed(batch_filepath, "<VALIDATE_RESULT>", str(validate).lower())

            ## - Walltime estimation:
            if mgcfd_unit_runtime_secs == 0.0:
                est_runtime_hours = 0
                est_runtime_minutes = 10
            else:
                est_runtime_secs = float(mgcfd_unit_runtime_secs*mg_cycles*mesh_multi) / math.sqrt(float(nt))
                est_runtime_secs *= 1.2 ## Allow for estimation error
                est_runtime_secs += 20  ## Add time for file load
                est_runtime_secs += 60
                est_runtime_secs = int(round(est_runtime_secs))
                est_runtime_hours = est_runtime_secs/60/60
                est_runtime_secs -= est_runtime_hours*60*60
                est_runtime_minutes = est_runtime_secs/60
                est_runtime_secs -= est_runtime_minutes*60
                if est_runtime_secs > 0:
                    est_runtime_minutes += 1
                    est_runtime_secs = 0
            py_sed(batch_filepath, "<HOURS>", str(est_runtime_hours).zfill(2))
            py_sed(batch_filepath, "<MINUTES>", str(est_runtime_minutes).zfill(2))

            ## Make batch script executable:
            os.chmod(batch_filepath, 0755)

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
            submit_all_file.write("\n\n")
            with open(submit_tmp_filepath, 'r') as f:
                for line in f:
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

    ## Create a little script for listing jobs that failed during execution:
    error_scan_filepath = os.path.join(jobs_dir, "list_errored_jobs.sh")
    error_scan_file = open(error_scan_filepath, "w")
    error_scan_file.write("#!/bin/bash\n")
    error_scan_file.write("find . -name pbs.stderr | while read F ; do wc -l \"$F\" ; done | grep -v \"^0 \"\n")
    error_scan_file.close()
    st = os.stat(error_scan_filepath)
    os.chmod(error_scan_filepath, st.st_mode | stat.S_IEXEC)
