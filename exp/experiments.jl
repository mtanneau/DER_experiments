import DataFrames
import CSV
using Statistics
using LinearAlgebra
using Random
using BenchmarkTools
using SparseArrays
using Printf

import MathProgBase
const MPB = MathProgBase

import CPLEX:CplexSolver
import Gurobi: GurobiSolver
import Mosek: MosekSolver

import Linda
import Tulip
BLAS.set_num_threads(1)

# Extend MPB 
import Mosek
import Gurobi
import CPLEX

MPB.getbarrieriter(m::Mosek.MosekMathProgSolverInterface.MosekLinearQuadraticModel) = Mosek.getintinf(m.task, Mosek.MSK_IINF_INTPNT_ITER)
MPB.getsimplexiter(m::Mosek.MosekMathProgSolverInterface.MosekLinearQuadraticModel) = Mosek.getintinf(m.task,Mosek.MSK_IINF_SIM_PRIMAL_ITER)+Mosek.getintinf(m.task,Mosek.MSK_IINF_SIM_DUAL_ITER)

MPB.getsimplexiter(m::Tulip.TulipMathProgModel) = 0

MPB.getsimplexiter(m::Gurobi.GurobiMathProgModel) = Gurobi.get_dblattr(m.inner, "IterCount")
MPB.getbarrieriter(m::Gurobi.GurobiMathProgModel) = Gurobi.get_intattr(m.inner, "BarIterCount")

MPB.getbarrieriter(m::CPLEX.CplexMathProgModel) = CPLEX.@cpx_ccall(getbaritcnt, Cint, (Ptr{Cvoid}, Ptr{Cvoid}), m.inner.env, m.inner.lp)
MPB.getsimplexiter(m::CPLEX.CplexMathProgModel) = CPLEX.@cpx_ccall(getitcnt, Cint, (Ptr{Cvoid}, Ptr{Cvoid}), m.inner.env, m.inner.lp)

# Demand Response module
include("../src/DemandResponse.jl")
# const DR = DemandResponse

"""
    generate_mp(n, T, netloadmin, netloadmax, cols, lpsolver)

Instanciate a Master Problem.

# Arguments
- `num_der`: The number of resources
- `T`: Length of the time-horizon (i.e. number of time-periods)
- `netloadmin`: Minimum total net load
- `netloadmax`: Maximum total net load
- `columns`: Initial set of columns
- `lpsolver`: LP solver for the RMP
"""
function generate_mp(
    num_der, T,
    netloadmin, netloadmax,
    columns,
    lpsolver
)
    # Generate RMP data
    A, obj, rowlb, rowub, collb, colub = generate_RMP_data(
        num_der, T, netloadmin, netloadmax,
        columns, lpsolver
    )

    # Instanciate RMP
    rmp = MPB.LinearQuadraticModel(lpsolver)
    MPB.loadproblem!(rmp, A, collb, colub, obj, rowlb, rowub, :Min)
    
    # Instanciate MP
    mp = Linda.LindaMaster(num_der, 2*T, vcat(netloadmax, netloadmin), 2*T, columns, rmp)
    
    return mp
end

"""
    generate_RMP_data(...)

Generate RMP data
"""
function generate_RMP_data(
    num_der, T, netloadmin, netloadmax,
    columns,
    lpsolver
)

    Alink = vcat(
        spzeros(num_der, 2*T),
        vcat(
            hcat(sparse(1.0I, T, T),  spzeros(T, T)), # <= u
            hcat(spzeros(T, T), sparse(-1.0I, T, T)) # >= l
        )
    )
    ncols = length(columns)
    A = hcat(
        Alink,
        vcat(
            hcat([sparsevec([c.idx_subproblem], [1.0], num_der) for c in columns]...),
            hcat([c.col for c in columns]...)
        )
    )
    
    obj = vcat(zeros(2*T), [c.cost for c in columns])
    
    rowlb = vcat(ones(num_der), netloadmax, netloadmin)
    rowub = vcat(ones(num_der), netloadmax, netloadmin)
    
    collb = zeros(2*T + ncols)
    colub = fill(Inf, 2*T + ncols)

    return A, obj, rowlb, rowub, collb, colub
end

function generate_RMP_data(
    num_der, T, netloadmin, netloadmax,
    columns,
    lpsolver::Tulip.TulipSolver
)
    ncols = length(columns)
    A = Tulip.TLPLinearAlgebra.UnitBlockAngular(
        Matrix(hcat([c.col for c in columns]...)),
        num_der,
        [c.idx_subproblem for c in columns],
        Matrix(vcat(
            hcat(sparse(1.0I, T, T),  spzeros(T, T)), # <= u
            hcat(spzeros(T, T), sparse(-1.0I, T, T)) # >= l
        ))
    )
    
    obj = vcat(zeros(2*T), [c.cost for c in columns])
    
    rowlb = vcat(ones(num_der), netloadmax, netloadmin)
    rowub = vcat(ones(num_der), netloadmax, netloadmin)
    
    collb = zeros(2*T + ncols)
    colub = fill(Inf, 2*T + ncols)

    return A, obj, rowlb, rowub, collb, colub
end

"""
    generate_resources()

# Arguments
- `num_der`: The number of DERs
- `T`: The number of time-steps
- `devices_rates`: Ownership rates of each device
- `load_norm`: Normalized reference load profile
- `PV_norm`: Normalized PV output profile
- `price`: Electricity price, in dollar per kW.h.
"""
function generate_resources(
    num_der, T,
    devices_rates::Dict{Symbol, Float64},
    load_norm::Vector,
    PV_norm::Vector,
    price::Vector;
    seed=0
)
    # Sanity checks
    T == length(load_norm) || throw(DimensionMismatch(
        "T=$T but load_norm has length $(length(load_norm))"
    ))
    T == length(PV_norm) || throw(DimensionMismatch(
        "T=$T but load_norm has length $(length(PV_norm))"
    ))
    T == length(price) || throw(DimensionMismatch(
        "T=$T but load_norm has length $(length(price))"
    ))

    Random.seed!(seed)
    resources = Vector{DR.DER.Resource}(undef, num_der)
    
    # Create a pool of `num_der` houses
    # Each household owns appliances according to `devices_rates`
    for i in 1:num_der

        appliances = DR.DER.Resource[]
        γ = 0.5 + rand()  # household's scaling factor

        # Fixed load
        l = DR.DER.FixedLoad(
            index=i,
            T=T,
            load= γ .* (load_norm .+ 0.05 .* rand(T))
        )
        push!(appliances, l)

        # Electric vehicle
        ndays = div(T, 24)
        if T==24
            ev = DR.DER.DeferrableLoad(
                index=i, T=T, dt=1.0,
                num_cycles=1,
                cycle_begin=[16],
                cycle_end=[24],
                cycle_energy_min=[10.0], cycle_energy_max=[10.0],
                pwr_min=[zeros(15); 1.1*ones(9)],
                pwr_max=[zeros(15); 7.7*ones(9)],
                binary_flag=true
            )
            push!(appliances, ev)
        end

        # PV
        ζ = rand(T)
        pv = DR.DER.CurtailableLoad(
            index=i,
            T=T,
            dt=1.0,
            load= γ .* PV_norm .* ζ,
            binaryFlag=true
        )
        push!(appliances, pv)

        # Battery
        bat = DR.DER.Battery(
            index=i,
            T=T,
            dt=1.0,
            soc_min=zeros(T), soc_max=10.0*ones(T), soc_init=0.0, selfdischargerate=1.0,
            charge_pwr_min=1.0, charge_pwr_max=6.0,
            discharge_pwr_min=1.0, discharge_pwr_max=6.0,
            charge_eff=1.0, discharge_eff=1.0
        )
        push!(appliances, bat)
        
        # Household
        h = DR.DER.House(
            index=i,
            T=T,
            dt=1.0,
            netload_min = zeros(T),
            netload_max = 10.0 .* ones(T),
            price=price,
            appliances=appliances
        )
        resources[i] = h
    end

    return resources
end