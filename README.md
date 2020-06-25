# DER experiments

This repository contains data and scripts for reproducing the computational experiments presented in Section 7 of the paper _Design and implementation of a modular interior-point solver for linear optimization_ (pre-print [here](https://arxiv.org/abs/2006.08814)).

## Installation instructions

1. Download and install Julia (version 1.3 or above)
2. Clone this repository
3. Download and install [Gurobi](https://www.gurobi.com)
4. Download and install Julia packages
    ```bash
    julia --project -e 'using Pkg; Pkg.instantiate();'
    ```

## Running experiments


To run the same set of experiments as in the paper, do the following:
1. Ensure the directories `exp/log` and `exp/rmp` exist and are empty.
2. Run individual jobs as follows:
    ```bash
    julia --project exp/colgen_grb.jl 24 1024
    ```
