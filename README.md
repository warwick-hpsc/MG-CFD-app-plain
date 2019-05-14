##########
# MG-CFD #
##########

An unstructured grid, multigrid (MG), finite-volume computational fluid dynamics (CFD) solver for inviscid flow. 
It has the goal of serving as a platform for evaluating emerging architectures, programming paradigms and algorithmic optimisations for this class of code.

This application is derived from Andrew Corrigan's CFD code as presented in the AIAA-2009-4001 paper, now included in the 'Rodinia' benchmark suite (http://rodinia.cs.virginia.edu/doku.php).

If you wish to cite this work then please refer to our (pending) journal publication:
Owenson A.M.B., Wright S.A., Bunt R.A., Ho Y.K., Street M.J., and Jarvis S.A. (2019), An Unstructured CFD Mini-Application for the Performance Prediction of a Production CFD Code, Concurrency Computat: Pract Exper., 2019

#########
# USAGE #
#########

1) Prepare a json file detailing run configuration. For example see ./run-inputs/simple.json
2) Generate run batch scripts using run-scripts/gen_job.py
3) Execute!
4) Collate together the output csv files using run-scripts/aggregate-output-data.py

###########
# CONTACT #
###########

Andrew Owenson: a.m.b.owenson@warwick.ac.uk
