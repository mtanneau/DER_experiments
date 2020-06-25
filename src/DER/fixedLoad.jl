"""
    FixedLoad

Fixed load, a.k.a. Uncontrollable load
"""
struct FixedLoad <: Resource
    index::Int           # Index of the resource
    num_timesteps::Int   # Number of time-steps in the battery's operation

    load::Vector{Float64}     # Fixed load at each time step

    function FixedLoad(;index::Int=0, T::Int=0, load::Vector{Float64}=Float64[])
        T == length(load) || throw(DimensionMismatch(
            "T=$T but load has size $(length(load))"
        ))
        return new(index, T, load)
    end
end

function add_resource_to_model!(
    m::MOI.ModelLike,
    h::House,
    l::FixedLoad,
    var2idx, con2idx
)
    T = l.num_timesteps

    # ==========================================
    #    I. Add local variables
    # ==========================================
    # Create variables in model
    pload = MOI.add_variables(m, T)
    # Set variable bounds
    # Indices of variable bounds are not tracked
    for t in 1:T
        MOI.add_constraint(m, MOI.SingleVariable(pload[t]), MOI.EqualTo(l.load[t]))
    end

    # Update list of variable indices
    for t in 1:T
        var2idx[(:fixed, l.index, :pnet, t)] = pload[t]  # net load
    end

    # ==========================================
    #    II. Update linking constraints
    # ==========================================
    for t in 1:T
        cidx = con2idx[(:house, h.index, :link, t)]
        MOI.modify(m, cidx, MOI.ScalarCoefficientChange(pload[t], -1.0))
    end

    # ==========================================
    #    III. Add local constraints
    # ==========================================
    # (None here)

    # ==========================================
    #    IV. Add sub-resources to current model
    # ==========================================
    # (None here)

    return nothing
end