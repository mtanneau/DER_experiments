# Numerical experiments

## Data

2016 hourly data from Ontario IESO, see files in `dat/` folder.

## Benchmark instances

* Time horizon: 24h, 48h, 72h, 96h
* Ownership ratesfor the devices

| Device            | Rate  |
|--                 | --:   |
| Fixed load        | 100%  |
| Dishwasher        | 65%   |
| Clothes washer    | 90%   |
| Clothes dryer     | 75%   |
| Heating           | 60%   |
| PV+EV+ES          | ~     |

* Proportion of renewables: `{0.00, 0.33, 0.66, 1.00}`
* Total aggregated load: in `[0.0, 6.0*R]`
* Number of resources: `{1024, 2048, 4096, 8192, 16384, 32768}`

## Solver for RMP

* IPM
    * Tulip
    * Mosek
    * CPLEX
    * Gurobi

* Simplex
    * S-CPLEX
    * S-Gurobi
    * S-CLP

## Parameters

* Partial pricing: 
    * stop pricing once `10% x R` have produced a violated cut
    * Full pricing is triggered if more than `25%` of sub-problems find no violated cut
* Reduced cost tolerance:
    * Primal / dual tolerance is `1e-6`
* Stopping criteria:
    * Duality gap < `1e-4`
    * Timeout at `2h` (only checked at the end of an iteration)



## Performance metrics

* CG: Number of outer iterations, number of columns, upper and lower bounds
* Computing times: RMP / SP / Total
* RMP: inner iterations
* SP: Number of calls to SP oracle (for information)

# Running the experiments

* For each instance
    * Generate instance and compute initial columns
    * For each solver
        * Create MP
        * Solve root node
        * Export results