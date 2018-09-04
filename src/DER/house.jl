"""
    House

Data stucture that represents a house and its appliances.

# Attributes
- `index`: Index of that house
- `num_timesteps`: Number of time-steps in the time horizon
- `dt`: Duration of one time-step, in hours
- `netload_min`: House's minimum net load, for each time-step
- `netload_max`:
- `price`: Electricity price, in $/kWh
- `appliances`: List of appliances in that house
"""
mutable struct House <: Resource
    index::Int

    num_timesteps::Int
    dt::Float64

    netload_min::Vector{Float64}
    netload_max::Vector{Float64}
    price::Vector{Float64}

    appliances::Vector{Resource}
end

function House(;
    index::Integer=0,
    T::Integer=0,
    dt::Float64=1.0,
    netload_min::Vector{Float64}=Vector{Float64}(),
    netload_max::Vector{Float64}=Vector{Float64}(),
    price::Vector{Float64}=Vector{Float64}(),
    appliances::Vector{Resource}=Vector{Resource}()
)

    # Dimension checks
    T == size(netload_min, 1) || throw(DimensionMismatch("Invalid netload_min"))
    T == size(netload_max, 1) || throw(DimensionMismatch("Invalid netload_max"))
    T == size(price, 1) || throw(DimensionMismatch("Invalid price"))
    0.0 < dt || throw(DomainError("dt must be positive"))

    House(
        index, T, dt,
        copy(netload_min), copy(netload_max), copy(price),
        appliances
    )
end

"""
    LindaOracleMIP(house, solver)

Construct a MIP model of the House, and casts it as a `LindaOracleMIP`.
The MIP model is built using JuMP, but only raw data is exported.
"""
function LindaOracleMIP(h::House, solver::MPB.AbstractMathProgSolver)
    
    T = h.num_timesteps
    m = JuMP.Model(solver=solver)  # Instanciate JuMP model

    #========================================
        Linking variables and constraints
    ========================================#
    netload = JuMP.@variable(
        m,
        [t=1:T],
        lowerbound=h.netload_min[t],
        upperbound=h.netload_max[t]
    )  # net load

    # minimize cost of purchased electricity
    JuMP.@objective(m, :Min, sum(h.price[t] * netload[t] for t=1:T))

    constr = JuMP.@constraint(
        m,
        [t=1:T],
        netload[t] == 0.0
    )  # net load balance
    
    
    #========================================
        Appliances
    ========================================#
    for a in h.appliances
        addmodel!(a, m, constr)
    end

    #========================================
        Export MPB model
    ========================================#
    JuMP.build(m)
    problem = JuMP.internalmodel(m)

    A_link = spzeros(2*h.num_timesteps, MPB.numvar(problem))
    for i in 1:h.num_timesteps
        A_link[i, i] = 1.0
        A_link[h.num_timesteps+i, i] = 1.0
    end

    return LindaOracleMIP(
        h.index,
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
    addmodel!(h::House, m::JuMP.Model, constr)

Add self to existing model.
"""
function addmodel!(h::House, m::JuMP.Model, constr)
    T = h.num_timesteps

    # Household's net load
    netload = [
        JuMP.@variable(
            m,
            objective=h.price[t],
            inconstraints=[constr[t]],
            coefficients=[-1.0],
            lowerbound=h.netload_min[t],
            upperbound=h.netload_max[t]
        )
        for t in 1:T
    ]

    # linking constraint at household level
    constr_ = JuMP.@constraint(
        m,
        [t=1:T],
        netload[t] == 0.0
    )

    # Recursively add appliances to model
    for a in h.appliances
        addmodel!(a, m ,constr_)
    end

    return
end