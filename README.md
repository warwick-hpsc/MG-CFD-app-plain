MG-CFD: An unstructured MG FVM CFD miniapp
==========================================

An unstructured grid, multigrid (MG), finite-volume computational fluid dynamics (CFD) solver for inviscid flow. 
It has the goal of serving as a platform for evaluating emerging architectures, programming paradigms and algorithmic optimisations for this class of code.

This application is derived from Andrew Corrigan's CFD code as presented in the AIAA-2009-4001 paper, now included in the 'Rodinia' benchmark suite (http://rodinia.cs.virginia.edu/doku.php).

Compiling and executing
==========================================

MG-CFD requires the [PAPI](https://icl.utk.edu/papi) library for collection of performance counter data.

### Instructions:

1) Prepare a json file detailing run configuration. See ./run-inputs/annotated.json for documentation on each option. 

2) Generate run batch scripts from the json file:

```Shell
     $ python ./run-scripts/gen_job.py --json path/to/config.json
```
     
3) The specified `jobs directory` will contain a subfolder for each run configuration, and a single `submit_all.sh` file. If a scheduler was specified in the .json file, then `submit_all.sh` will compile locally then submit each job to the scheduler for execution. If local execution was requested in the json file, then `submit_all.sh` will compile and execute locally. 

4) Each run will output .csv files containing performance data, as well as the compiler-generated object files for performance-critical loops. These can be collated together using `aggregate-output-data.py`. [assembly-loop-extractor](https://github.com/warwick-hpsc/assembly-loop-extractor) is required to analyse the object files.

```Shell
     $ python ./run-scripts/aggregate-output-data.py \
              --assembly-dirpath path/to/git-clone-of/assembly-loop-extractor \
              --output-dirpath path/to/desired-folder/for/collated-csv-files \
              --data-dirpaths path/to/job-group-1-output [path/to/job-group-2-output ...]
```

Datasets
==========================================

A release is provided that includes two meshes. The first is the `fvcorr.domn.097K` originally bundled with the original CFD code, enabling numerical validation between that and CFD. 

The second mesh is of the [Onera M6 wing](https://www.grc.nasa.gov/WWW/wind/valid/m6wing/m6wing.html). It consists of 300K nodes (930K edges), and three additional multigrid meshes with respective node counts of 165K, 111K, and 81K.

Authorship
==========================================

Andrew Owenson: a.m.b.owenson@warwick.ac.uk

If you wish to cite this work then please refer to our (pending) journal publication:

* Owenson A.M.B., Wright S.A., Bunt R.A., Ho Y.K., Street M.J., and Jarvis S.A. (2019), An Unstructured CFD Mini-Application for the Performance Prediction of a Production CFD Code, Concurrency Computat: Pract Exper., 2019

Release
==========================================

MG-CFD is released under a MIT license.
