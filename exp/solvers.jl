import MathProgBase
const MPB = MathProgBase

# Import solvers
import Linda
import Tulip
import Mosek
import Gurobi
import CPLEX

import CPLEX:CplexSolver
import Gurobi: GurobiSolver
import Mosek: MosekSolver
import Tulip.TulipSolver


# Extend MathProgBase interface
MPB.getbarrieriter(m::CPLEX.CplexMathProgModel) = CPLEX.@cpx_ccall(getbaritcnt, Cint, (Ptr{Cvoid}, Ptr{Cvoid}), m.inner.env, m.inner.lp)
MPB.getsimplexiter(m::CPLEX.CplexMathProgModel) = CPLEX.@cpx_ccall(getitcnt, Cint, (Ptr{Cvoid}, Ptr{Cvoid}), m.inner.env, m.inner.lp)

MPB.getsimplexiter(m::Gurobi.GurobiMathProgModel) = Gurobi.get_dblattr(m.inner, "IterCount")
MPB.getbarrieriter(m::Gurobi.GurobiMathProgModel) = Gurobi.get_intattr(m.inner, "BarIterCount")

MPB.getbarrieriter(m::Mosek.MosekMathProgSolverInterface.MosekLinearQuadraticModel) = Mosek.getintinf(m.task, Mosek.MSK_IINF_INTPNT_ITER)
MPB.getsimplexiter(m::Mosek.MosekMathProgSolverInterface.MosekLinearQuadraticModel) = Mosek.getintinf(m.task,Mosek.MSK_IINF_SIM_PRIMAL_ITER)+Mosek.getintinf(m.task,Mosek.MSK_IINF_SIM_DUAL_ITER)

MPB.getsimplexiter(m::Tulip.TulipMathProgModel) = 0

# Solver instanciation

solver(s::Symbol) = solver(Val(s))
solver(::Val{:TLP}) = generate_lpsolver(:TLP)
solver(::Val{:MSK}) = generate_lpsolver(:MSK, lpmethod=4, crossover=false)
solver(::Val{:GRB}) = generate_lpsolver(:GRB, lpmethod=4, crossover=false)
solver(::Val{:CPX}) = generate_lpsolver(:CPX, lpmethod=4, crossover=false)
solver(::Val{:CPX_}) = generate_lpsolver(:CPX_, lpmethod=4, crossover=false)
solver(::Val{:SGRB}) = generate_lpsolver(:GRB, lpmethod=1, crossover=true)
solver(::Val{:SCPX}) = generate_lpsolver(:CPX, lpmethod=1, crossover=true)


"""
    generate_lpsolver()

Instanciate an LP solver for the RMP.

# Arguments
- `name`: Name of the solver
- `verbose`: verbosity level. 0 means no output.
- `lpmethod`: indicates which algorithm to use.
    0 is default, 1 is primal simplex, 2 is dual simplex, 3 is simplex, 4 is barrier
- `nthreads`: number of threads available to the solver
- `crossover`: whether to perform crossover.
    0 is default, 1 is yes, 2 is no
- `presolve`: whether to perform presolve
    0 is no, 1 is yes
"""
generate_lpsolver(
    name::Symbol;
    verbose::Int=0,
    lpmethod::Int=0,
    nthreads=1,
    crossover::Bool=false,
    presolve::Int=0
) = generate_lpsolver(Val(name), verbose, lpmethod, nthreads, crossover, presolve)


function generate_lpsolver(::Val{:CPX},
    _verbose, _lpmethod, _nthreads, _crossover, _presolve
)

    # crossover parameter
    if _lpmethod == 4 
        crossover = _crossover ? 1 : 2
    else
        crossover = 0
    end

    CplexSolver(
        CPX_PARAM_SCRIND=_verbose,
        CPX_PARAM_PREIND=_presolve,
        CPX_PARAM_THREADS=_nthreads,
        CPX_PARAM_LPMETHOD=_lpmethod, # 0:default, 1:PS, 2:DS, 4:IPM, 6:Concurrent
        CPX_PARAM_SOLUTIONTYPE=crossover  # 0: default, 1: basic, 2:no crossover,
    )

end

function generate_lpsolver(::Val{:GRB},
    _verbose, _lpmethod, _nthreads, _crossover, _presolve
)

    # crossover parameter
    crossover = _crossover ? 1 : 0

    if _lpmethod == 1 || _lpmethod == 2
        lpmethod = _lpmethod - 1  # simplex
    elseif _lpmethod == 4
        lpmethod = 2  # barrier
    else
        lpmethod = -1  # default setting
    end

    GurobiSolver(
        OutputFlag=_verbose,
        Method=lpmethod,
        Presolve=_presolve,
        Threads=_nthreads,
        Crossover=crossover,
        InfUnbdInfo=1  # to get infeasible rays if needed
    )

end

function generate_lpsolver(::Val{:MSK},
    _verbose, _lpmethod, _nthreads, _crossover, _presolve
)
    if _lpmethod == 1
        optimizer = 6  # primal simplex
    elseif _lpmethod == 2
        optimizer = 1  # dual simplex
    elseif _lpmethod == 3
        optimizer = 3  # any simplex
    elseif _lpmethod == 4
        optimizer = 4  # interior-point
    else
        optimizer = 2  # default
    end
    MosekSolver(
        MSK_IPAR_LOG=_verbose,
        MSK_IPAR_NUM_THREADS=_nthreads,
        MSK_IPAR_OPTIMIZER=optimizer,
        MSK_IPAR_PRESOLVE_USE=_presolve,
        MSK_IPAR_INTPNT_BASIS=(1*_crossover),  # crossover
        # MSK_IPAR_INTPNT_MAX_NUM_COR=-1  # default
    )

end

generate_lpsolver(::Val{:TLP},
    _verbose, _lpmethod, _nthreads, _crossover, _presolve
) = TulipSolver(verbose = _verbose, barrier_iter_max=200)