mutable struct House <: Resource
    index::Int

    num_timesteps::Int
    dt::Float64

    netload_min::Vector{Float64}
    netload_max::Vector{Float64}
    price::Vector{Float64}

    appliances::Vector{Resource}
end

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