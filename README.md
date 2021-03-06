MG-CFD: An unstructured MG FVM CFD miniapp
==========================================

A 3D unstructured multigrid, finite-volume computational fluid dynamics (CFD) mini-app for inviscid-flow. 
It has the goal of serving as a platform for evaluating emerging architectures, programming paradigms and algorithmic optimisations for this class of code. 

This application is derived from Andrew Corrigan's `CFD Solver` code as presented in the AIAA-2009-4001 paper, now included in the `Rodinia` benchmark suite (http://rodinia.cs.virginia.edu/doku.php).

Compiling and executing
==========================================

1) Clone this repository

2) Pick a target compiler: gnu, intel, clang, or cray

3) Compile: `COMPILER=intel make`

4) Navigate to input data folder and quick run: `path/to/bin/euler3d*.b -i input.dat`

MG-CFD has more command-line arguments to ease file/directory interaction, and control execution parameters. View the help page for more information:

```Shell
     $ ./path/to/euler3d*.b --help
```

### Generating batch submission scripts:

MG-CFD comes with templates and logic for generating submission scripts to one of several job schedulers, that can iterate over combinations of MG-CFD parameters. To generate them follow these instructions:

1) Prepare a json file detailing run configuration. See [run-inputs/annotated.json](https://github.com/warwick-hpsc/MG-CFD-app-plain/blob/master/run-inputs/annotated.json) for documentation on each option. 

2) Generate run batch scripts from the json file:

```Shell
     $ python ./run-scripts/gen_job.py --json path/to/config.json
```
     
3) The specified `jobs directory` will contain a subfolder for each run configuration, and a single `submit_all.sh` file. If a scheduler was specified in the .json file, then `submit_all.sh` will compile locally then submit each job to the scheduler for execution. If local execution was requested in the json file, then `submit_all.sh` will compile and execute locally. 

4) Each run will output CSV files containing performance data, as well as the compiler-generated object files for performance-critical loops. These can be collated together using `aggregate-output-data.py`:

```Shell
     $ python ./run-scripts/aggregate-output-data.py \
              --output-dirpath path/to/desired-folder/for/collated-csv-files \
              --data-dirpaths path/to/job-group-1-output [path/to/job-group-2-output ...]
```

5) To also analyse the object files, [assembly-loop-extractor](https://github.com/warwick-hpsc/assembly-loop-extractor) is required:

```Shell
     $ python ./run-scripts/aggregate-output-data.py \
              --assembly-dirpath path/to/git-clone-of/assembly-loop-extractor \
              --output-dirpath path/to/desired-folder/for/collated-csv-files \
              --data-dirpaths path/to/job-group-1-output [path/to/job-group-2-output ...]
```

Performance assessment
==========================================

One role of MG-CFD is to evaluate a single compute node in regards to performance of this class of code. This involves measuring the throughputs of floating-point arithmetic, loads & stores, and main memory performance. To collect this data, two script generation templates have been provided with MG-CFD, named `assess-compute.json` and `assess-memory.json`. 

To assess compute performance, MG-CFD requires the [PAPI](https://icl.utk.edu/papi) library for collection of performance counter data (counts of instructions and clock cycles).

Once collected and aggregated, this data can be passed into the [MG-CFD performance model](https://github.com/warwick-hpsc/MG-CFD-performance-model) to estimate the throughput rates.

Datasets
==========================================

A [release](https://github.com/warwick-hpsc/MG-CFD-app-plain/releases) is provided that includes two meshes. The first is the `fvcorr.domn.097K` originally bundled with the original `CFD Solver` code, enabling numerical validation between it and `MG-CFD`. 

The second mesh is of the [Onera M6 wing](https://www.grc.nasa.gov/WWW/wind/valid/m6wing/m6wing.html). It consists of 300K nodes (930K edges), and three additional multigrid meshes with respective node counts of 165K, 111K, and 81K.

Authorship
==========================================

Andrew Owenson: a.m.b.owenson@warwick.ac.uk

For more information on the design of MG-CFD and performance model, please refer to our publication: https://onlinelibrary.wiley.com/doi/10.1002/cpe.5443

If you wish to cite this work then please use the following:

* Owenson A.M.B., Wright S.A., Bunt R.A., Ho Y.K., Street M.J., and Jarvis S.A. (2019), An Unstructured CFD Mini-Application for the Performance Prediction of a Production CFD Code, Concurrency Computat: Pract Exper., 2019

