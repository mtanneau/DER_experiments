"""
    FixedLoad

Fixed load, a.k.a. Uncontrollable load
"""
mutable struct FixedLoad <: Resource
    index::Int           # Index of the resource
    num_timesteps::Int   # Number of time-steps in the battery's operation

    load::Vector{Float64}     # Fixed load at each time step
end

function FixedLoad(;
    index::Integer=0,
    T::Integer=0,
    dt::Float64=1.0,
    load::AbstractVector{T1}=[0.0]
) where{T1<:Real}

    # Dimension checks
    T == size(load, 1) || throw(DimensionMismatch("Invalid soc_min"))
    0.0 < dt || throw(DomainError("dt must be positive"))

    FixedLoad(index, T, load)
end

"""
    LindaOracleMIP(l::FixedLoad, solver::AbstractMathProgSolver)

Construct a MIP oracle that schedules the Battery.
"""
function LindaOracleMIP(l::FixedLoad, solver::MPB.AbstractMathProgSolver)
    T = l.num_timesteps
    # Construct model
    m = JuMP.Model(solver=solver)

    # add variables
    JuMP.@variable(m, p[t=1:T]==l.load[t])  # net power
    
    # set objective
    # (minimize total cost of purchased energy)
    JuMP.@objective(m, :Min, 0.0)

    JuMP.build(m)
    problem = JuMP.internalmodel(m)

    A_link = spzeros(2*bat.num_timesteps, MPB.numvar(problem))
    for i in 1:bat.num_timesteps
        A_link[i, i] = 1.0
        A_link[bat.num_timesteps+i, i] = 1.0
    end

    return LindaOracleMIP(
        bat.index,
        MPB.getobj(problem),
        A_link,
        MPB.getconstrmatrix(problem),
        MPB.getconstrLB(problem),
        MPB.getconstrUB(problem),
        MPB.getvartype(problem),
        MPB.getvarLB(problem),
        MPB.getvarUB(problem),
        solver
    )
end

"""
    addmodel!(l::FixedLoad, m::JuMP.Model, constr)

Add battery to existing model.
"""
function addmodel!(l::FixedLoad, m::JuMP.Model, constr)
    T = l.num_timesteps

    #========================================
        Local variables
    ========================================#
    # Net power, appears in linking constraints with coefficient `-1.0`
    p = [
        JuMP.@variable(
            m,
            objective=0.0,
            inconstraints=[constr[t]],
            coefficients=[-1.0],
            lowerbound=l.load[t],
            upperbound=l.load[t]
        )
        for t in 1:T
    ]

    return
end