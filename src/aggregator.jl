import JuMP

import MathProgBase
const MPB = MathProgBase

"""
    Aggregator
"""    
mutable struct Aggregator
    num_timesteps::Int   # Number of time steps

    netload_min::Vector{Float64}    # Minimum net load
    netload_max::Vector{Float64}    # Maximum net load

    price_mkt::Vector{Float64}  # Market price for electricity
    resources::Vector{DER.Resource} # List of DERs
end

function build_centralized_model(a::Aggregator, solver::MPB.AbstractMathProgSolver)
    T = a.num_timesteps

    # Instanciate model
    m = JuMP.Model(solver=solver)

    # Linking variables
    netload = JuMP.@variable(
        m,
        [t=1:T],
        lowerbound=a.netload_min[t],
        upperbound=a.netload_max[t],
    )

    # Minimize cost of purchased electricity
    JuMP.@objective(
        m,
        :Min,
        sum(a.price_mkt[t] * netload[t] for t=1:T)
    )

    # Linking constraints
    constr_link = JuMP.@constraint(
        m,
        [t=1:T],
        netload[t] == 0.0
    )  # net load balance

    # Add houses to model
    for r in a.resources
        DER.addmodel!(r, m, constr_link)
    end

    JuMP.build(m)

    return JuMP.internalmodel(m)

end