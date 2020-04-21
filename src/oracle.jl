import MathOptInterface
const MOI = MathOptInterface

import Linda:
    Column,
    Oracle, update!, optimize!, get_columns, get_dual_bound, get_objective_value
const House = DER.House

mutable struct HHOracle <: Oracle

    h::House
    model::MOI.ModelLike

    var2idx::Dict{Tuple{Symbol, Int, Symbol, Int}, MOI.VariableIndex}    # Variable indices in the model
    con2idx::Dict{Tuple{Symbol, Int, Symbol, Int}, MOI.ConstraintIndex}  # Constraint indices in the model

    # Shadow prices
    farkas::Bool
    π::Vector{Float64}
    σ::Float64
end

function update!(o::HHOracle, farkas::Bool, π::Vector{Float64}, σ::Float64)

    # Sanity checks
    T = o.h.num_timesteps
    length(π) == 2*T || throw(DimensionMismatch("π has length $(length(π)) but should have length 2xT=$(2*T)."))
    
    # Update shadow prices
    o.farkas = farkas
    o.π .= π[1:T] .+ π[T+1:end]
    o.σ  = σ

    # Update objective coefficients
    pnet = [
        o.var2idx[(:house, o.h.index, :pnet, t)]
        for t in 1:T
    ]
    for t in 1:T
        MOI.modify(o.model,
            MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
            MOI.ScalarCoefficientChange(pnet[t], (!farkas) * o.h.price[t] - o.π[t])
        )
    end

    return nothing
end

function optimize!(o::HHOracle)
    MOI.optimize!(o.model)
    return nothing
end

function get_columns(o::HHOracle)

    T = o.h.num_timesteps

    # TODO: check solution status
    pnet = [
        MOI.get(o.model, MOI.VariablePrimal(),
            o.var2idx[(:house, o.h.index, :pnet, t)]
        )
        for t in 1:T
    ]

    # Compute real cost
    cost = dot(o.h.price, pnet)

    # Get reduced cost
    rc = MOI.get(o.model, MOI.ObjectiveValue()) - o.σ

    # Create column
    col = Column(
        cost,
        o.h.index,
        collect(1:2*T), [pnet; pnet],
        true
    )

    return [(col, rc)]
end

function get_dual_bound(o::HHOracle)
    if o.farkas
        return -Inf
    else
        return MOI.get(o.model, MOI.ObjectiveValue()) - o.σ
    end
end

# Create an oracle from a house
function build_oracle(m::MOI.ModelLike, h::House)
    # Build pricing problem
    _, var2idx, con2idx = DR.DER.build_model!(h, m)

    # Setup oracle
    T = h.num_timesteps
    return HHOracle(h, m, var2idx, con2idx, false, zeros(T), 0.0)
end