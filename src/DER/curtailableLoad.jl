"""
    CurtailableLoad

Curtailable load
"""
struct CurtailableLoad <: Resource
    index::Int           # Index of the resource
    num_timesteps::Int   # Number of time-steps in the battery's operation

    load::Vector{Float64}     # Load at each time-step, if device is on
    binflag::Bool  # Whether on-off or continuous curtailment

    function CurtailableLoad(;
        index::Integer=0,
        T::Integer=0,
        load::Vector{Float64}=Float64[],
        binflag::Bool=true
    )
        # Dimension checks
        T == size(load, 1) || throw(DimensionMismatch("Invalid load dimension"))
        return new(index, T, load, binflag)
    end

end

function add_resource_to_model!(
    m::MOI.ModelLike,
    h::House,
    l::CurtailableLoad,
    var2idx, con2idx
)
    T = l.num_timesteps

    # ==========================================
    #    I. Add local variables
    # ==========================================
    pnet = MOI.add_variables(m, T)  # Net load
    curt = MOI.add_variables(m, T)  # Curtailment indicator
    # Update list of variable indices
    for t in 1:T
        var2idx[(:curt, l.index, :pnet, t)] = pnet[t]
        var2idx[(:curt, l.index, :uind, t)] = curt[t]
    end

    for t in 1:T
        if l.binflag
            # Add binary constraint
            MOI.add_constraint(m,
                MOI.SingleVariable(curt[t]),
                MOI.ZeroOne()
            )
        else
            # Just add bounds
            MOI.add_constraint(m,
                MOI.SingleVariable(curt[t]),
                MOI.Interval(0.0, 1.0)
            )
        end
    end

    # ==========================================
    #    II. Update linking constraints
    # ==========================================
    for t in 1:T
        cidx = con2idx[(:house, h.index, :link, t)]
        MOI.modify(m, cidx, MOI.ScalarCoefficientChange(pnet[t], -1.0))
    end

    # ==========================================
    #    III. Add local constraints
    # ==========================================
    # Curtailement
    #   p[t] - u[t] * P[t] = 0
    for t in 1:T
        cidx = MOI.add_constraint(m,
            MOI.ScalarAffineFunction{Float64}(
                [
                    MOI.ScalarAffineTerm(1.0, pnet[t]),
                    MOI.ScalarAffineTerm(-l.load[t], curt[t])
                ],
                0.0
            ),
            MOI.EqualTo(0.0)
        )
        con2idx[(:curt, l.index, :curtail, t)] = cidx
    end
    
    # ==========================================
    #    IV. Add sub-resources to current model
    # ==========================================
    # (None here)

    return nothing
end