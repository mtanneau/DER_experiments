# DER experiments

This repository contains data and scripts for reproducing the computational experiments of 

## Installation instructions

### Install Julia

### Install solvers

1. Commercial solvers require a license and must be installed separately.
    * [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio)
    * [Gurobi](https://www.gurobi.com)
    * [Mosek](https://www.mosek.com)

2.Tulip will be downloaded automatically when installing Julia packages.

### Clone this repository

```bash
git clone https://github.com/mtanneau/DER_experiments.git
```

### Install Julia packages

All package dependencies are specified in the `Project.toml` and `Manifest.toml` files located at the root of this directory.

To download and install those packages, run Julia from the root of this directory, and run the following commands
```julia
> using Pkg
> Pkg.activate(".")
> Pkg.instantiate()
```
or, from the terminal:
```bash
julia --project=. -e 'using Pkg; Pkg.instantiate();' 
```

Additional steps may be required for building the Julia interface of commercial solvers. See installation instructions for each solver's API.

## Running experiments

Each individual experiment is performed by running the `benchmark.jl` script.
The correspoding syntax is:
```bash
julia --project=. scripts/benchmark.jl -T 24 -R 128 -s 0 --export
```
The meaning of each parameter is explained in the code.


To run the same set of experiments as in the paper, do the following:

1. Generate the list of jobs by running
```bash
julia --project=. scripts/exp_preprocessing.jl
```
This will create a `jobs.txt` file that contains the list of all individual jobs.

2. Create a `res/` and `out` directories 
```bash
mkdir scripts/res
mkdir scripts/out
```
If such directories already exist, make sure to empty them before running the experiments.

3. Run each job, e.g. using GNU parallel.
Note that large experiments may take several hours to run (up to 10,000s for each of the solver).

## Post-processing

For each individual experiment (i.e., each tuple `(T, R, s)`), the corresponding results file will be located in `scripts/res`.
Each file has the form `cg_<T>_<R>_0.33_<s>_<solver>.csv`, e.g.
```
cg_24_1024_0.33_5_CPX.csv
```