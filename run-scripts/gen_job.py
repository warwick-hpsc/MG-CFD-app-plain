import os, shutil, stat
import json, argparse, re
import itertools
import math

script_dirpath = os.path.dirname(os.path.realpath(__file__))
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
js_to_submit_cmd["pbs"] = "qsub"

defaults = {}
defaults["compiler"] = "intel"
defaults["insn sets"] = [ "Host"]
defaults["base flags"] = "-DTIME"
defaults["perm flags"] = [""]
defaults["num threads"] = [1]
defaults["num repeats"] = 1
defaults["mg cycles"] = 50
defaults["min mesh multi"] = 1
## Approximate runtime of a single MG cycle of LA_cascade mesh:
defaults["unit walltime"] = 0.2
defaults["budget code"] = "NotSpecified"

def get_key_value(profile, cat, key):
    if not cat in profile.keys():
        raise Exception("Cat '{0}' not in json.".format(cat))

    if key in profile[cat]:
        return profile[cat][key]
    else:
        if key in defaults:
            return defaults[key]
        else:
            raise Exception("Mandatory key '{0}' not present in cat '{1}' of json".format(key, cat))

def py_sed(filepath, from_rgx, to_rgx):
    with open(filepath, "r") as f:
        lines = f.readlines()
    with open(filepath, "w") as f:
        for line in lines:
            if isinstance(to_rgx, str):
                f.write(re.sub(from_rgx,     to_rgx,  line))
            else:
                f.write(re.sub(from_rgx, str(to_rgx), line))

def prune_perm_flags_permutations(perm_flags_permutations):
    perm_flags_permutations_pruned = []
    for p in perm_flags_permutations:
        flux_cripple_option_present = False
        other_flux_options_present = False
        for pi in p:
            if "FLUX_CRIPPLE" in pi:
                flux_cripple_option_present = True
            elif "FLUX_" in pi:
                other_flux_options_present = True
        if flux_cripple_option_present and other_flux_options_present:
            ## Discard this permutation, as the FLUX_CRIPPLE option should be
            ## applied without other flux options.
            pass
        else:
            perm_flags_permutations_pruned.append(p)
    return perm_flags_permutations_pruned

def delete_folder_contents(dirpath):
    print("Deleting contents of folder: " + dirpath)
    for f in os.listdir(dirpath):
        fp = os.path.join(dirpath, f)
        if os.path.isdir(fp):
            shutil.rmtree(fp)
        else:
            os.remove(fp)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--json', required=True)
    args = parser.parse_args()
    with open(args.json, 'r') as f:
        profile = json.load(f)

    jobs_dir = profile["setup"]["jobs dir"]
    if jobs_dir[0] != '/':
        jobs_dir = os.path.join(app_dirpath, jobs_dir)
    if not os.path.isdir(jobs_dir):
        os.mkdir(jobs_dir)
    else:
        delete_folder_contents(jobs_dir)

    data_dirpath = profile["run"]["data dirpath"]
    if data_dirpath[0] != '/':
        data_dirpath = os.path.join(app_dirpath, data_dirpath)

    submit_all_filepath = os.path.join(jobs_dir, "submit_all.sh")
    submit_all_file = open(submit_all_filepath, "w")
    submit_all_file.write("#!/bin/bash\n\n")

    js = profile["setup"]["job scheduler"]
    js_filename = js_to_filename[js]
    submit_all_file.write("# {0}:\n".format(js))
    submit_all_file.write("submit_cmd={0}\n\n".format(js_to_submit_cmd[js]))

    job_queue = get_key_value(profile, "setup", "partition")
    budget_code = get_key_value(profile, "setup", "budget code")

    with open(os.path.join(jobs_dir, "papi.conf"), "w") as f:
        f.write("PAPI_TOT_INS\n")
        f.write("PAPI_TOT_CYC\n")

    isas = get_key_value(profile, "compile", "insn sets")
    threads = get_key_value(profile, "run", "num threads")
    num_repeats = get_key_value(profile, "run", "num repeats")
    mg_cycles = get_key_value(profile, "run", "mg cycles")
    mgcfd_unit_runtime_secs = get_key_value(profile, "run", "unit walltime")
    compiler = get_key_value(profile, "compile", "compiler")
    base_flags = get_key_value(profile, "compile", "base flags")
    perm_flags = get_key_value(profile, "compile", "perm flags")
    perm_flags_permutations = prune_perm_flags_permutations(itertools.product(*perm_flags))
    if len(perm_flags_permutations) == 0:
        perm_flags_permutations = [""]

    min_mesh_multi = get_key_value(profile, "run", "min mesh multi")

    num_jobs = len(perm_flags_permutations) * len(isas) * len(threads) * num_repeats
    submit_all_file.write("num_jobs={0}\n\n".format(num_jobs))

    n = 0
    for p in perm_flags_permutations:
        flags = ""
        for pi in p:
            if pi != "":
                if flags != "":
                    flags += " "
                flags += pi
        for isa in isas:
            for nt in threads:
                for repeat in range(num_repeats):
                    n += 1
                    job_id = str(n).zfill(3)

                    print("Creating job {0}/{1}".format(n, num_jobs))

                    job_dir = os.path.join(jobs_dir, job_id)
                    if not os.path.isdir(job_dir):
                        os.mkdir(job_dir)

                    ## Prepare MG-CFD execution
                    build_flags = base_flags + " " + flags

                    dest_filepath = os.path.join(job_dir, "papi.conf")
                    if os.path.isfile(dest_filepath):
                        os.remove(dest_filepath)
                    os.symlink(os.path.join(jobs_dir, "papi.conf"), dest_filepath)

                    job_run_filepath = os.path.join(job_dir, "run-mgcfd.sh")
                    shutil.copyfile(os.path.join(template_dirpath, "run-mgcfd.sh"), job_run_filepath)
                    py_sed(job_run_filepath, "<RUN_OUTDIR>", job_dir)
                    py_sed(job_run_filepath, "<APP_DIRPATH>", app_dirpath)
                    py_sed(job_run_filepath, "<DATA_DIRPATH>", data_dirpath)
                    py_sed(job_run_filepath, "<ISA>", isa)
                    py_sed(job_run_filepath, "<BUILD_FLAGS>", build_flags)
                    py_sed(job_run_filepath, "<COMPILER>", compiler)
                    py_sed(job_run_filepath, "<NUM_THREADS>", nt)
                    py_sed(job_run_filepath, "<MG_CYCLES>", mg_cycles)

                    while (mesh_multi % nt) > 0:
                        mesh_multi += 1
                    py_sed(job_run_filepath, "<MESH_MULTI>", mesh_multi)

                    ## Prepare job scheduling:
                    if js != "":
                        js_filepath = os.path.join(job_dir, js_filename)
                        shutil.copyfile(os.path.join(template_dirpath, js_filename), js_filepath)
                        py_sed(js_filepath, "<PARTITION>", job_queue)
                        py_sed(js_filepath, "<RUN_DIR>", job_dir)
                        py_sed(js_filepath, "<BUDGET CODE>", budget_code)

                        est_runtime_secs = float(mgcfd_unit_runtime_secs*mg_cycles*mesh_multi) / math.sqrt(float(nt))
                        est_runtime_secs = 1.2*est_runtime_secs + 10.0 ## Add a small buffer
                        est_runtime_secs = int(round(est_runtime_secs))
                        est_runtime_hours = est_runtime_secs/60/60
                        est_runtime_secs -= est_runtime_hours*60*60
                        est_runtime_minutes = est_runtime_secs/60
                        est_runtime_secs -= est_runtime_minutes*60
                        if est_runtime_secs > 0:
                            est_runtime_minutes += 1
                            est_runtime_secs = 0
                        py_sed(js_filepath, "<HOURS>", str(est_runtime_hours).zfill(2))
                        py_sed(js_filepath, "<MINUTES>", str(est_runtime_minutes).zfill(2))

                    ## Combine into a batch submission script:
                    if js == "":
                        batch_filename = "run.sh"
                    else:
                        batch_filename = js+".batch"
                    batch_filepath = os.path.join(job_dir, batch_filename)
                    with open(batch_filepath, "w") as f_out:
                        if js != "":
                            with open(js_filepath, "r") as f_in:
                                for line in f_in.readlines():
                                    f_out.write(line)
                            os.remove(js_filepath)
                            f_out.write("\n")
                        with open(job_run_filepath, "r") as f_in:
                            for line in f_in.readlines():
                                f_out.write(line)
                        os.remove(job_run_filepath)

                    ## Make batch script executable:
                    os.chmod(batch_filepath, 0755)

                    ## Add an entry to submit_all.sh:
                    submit_all_file.write("\n\n")
                    submit_all_file.write('if [ ! -f "{0}"/Times.csv ]; then\n'.format(job_dir))
                    submit_all_file.write("  basedir=`pwd`\n")
                    submit_all_file.write("  cd {0}\n".format(job_dir))
                    if js == "":
                        ## Without a job scheduler to log STDOUT, need to do this manually:
                        submit_all_file.write('  echo "Executing job {0}/{1}"\n'.format(n, num_jobs))
                        submit_all_file.write('  eval "$submit_cmd" ./{0} > submit.log\n'.format(batch_filename))
                    else:
                        submit_all_file.write('  eval "$submit_cmd" ./{0}\n'.format(batch_filename))
                    submit_all_file.write('  cd "$basedir"\n')
                    submit_all_file.write('fi')

    submit_all_file.close()
    st = os.stat(submit_all_filepath)
    os.chmod(submit_all_filepath, st.st_mode | stat.S_IEXEC)