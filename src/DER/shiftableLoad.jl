"""
    ShiftableLoad

ShiftableLoad Load
"""
struct ShiftableLoad <: Resource
    index::Int  # Unique index of device
    num_timesteps::Int  # Number of time-steps in horizon

    cycle_begin_min::Int  # Beginning of each cycle
    cycle_begin_max::Int    # End of each cycle
    cycle_length::Int  # Length of a cycle
    cycle_load::Vector{Float64}  # Load level at each time of the cycle

    function ShiftableLoad(;
        index::Integer=0,
        T::Integer=0,
        cycle_begin_min::Int=1,
        cycle_begin_max::Int=0,
        cycle_length::Int=0,
        cycle_load::Vector{Float64}=Float64[]
    )
        # Dimension checks
        T >= cycle_length || throw(error(
            "T=$T but cycle has length $cycle_length"
        ))
        cycle_begin_min >= 1 || throw(error(
            "Earliest start time is too early: $cycle_begin_min"
        ))
        cycle_begin_max <= T - cycle_length + 1 || throw(error(
            "Latest start time is too late: $cycle_begin_max"
        ))
        cycle_length == length(cycle_load) || throw(DimensionMismatch(
            "L=$cycle_length but load has length $(length(cycle_load))"
        ))

        return new(
            index, T,
            cycle_begin_min, cycle_begin_max,
            cycle_length, cycle_load
        )
    end
end

function add_resource_to_model!(
    m::MOI.ModelLike,
    h::House,
    l::ShiftableLoad,
    var2idx, con2idx
)
    T = l.num_timesteps
    T0 = l.cycle_begin_min
    T1 = l.cycle_begin_max

    # ==========================================
    #    I. Add local variables
    # ==========================================
    pnet = MOI.add_variables(m, T)  # Net load
    L = T1 - T0 + 1
    ustart = MOI.add_variables(m, L)  # Start indicator

    # Variable bounds
    for k in 1:L
        MOI.add_constraint(m, MOI.SingleVariable(ustart[k]), MOI.ZeroOne())
    end

    # Update variable indices
    for t in 1:T
        var2idx[(:shift, l.index, :pnet, t)] = pnet[t]
    end
    for k in 1:L
        var2idx[(:shift, l.index, :pnet, k + T0 - 1)] = ustart[k]
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
    # Net load
    for t in 1:T
        MOI.add_constraint(m, 
            MOI.ScalarAffineFunction([
                [
                    MOI.ScalarAffineTerm(1.0, pnet[t])
                ]; [
                    MOI.ScalarAffineTerm(-l.cycle_load[t - T0 - k + 2], ustart[k])
                    for k in 1:L
                    if (T0 <= t - k + 1 <= T1) && (1 <= t - T0 - k + 2 <= l.cycle_length)
                ]
            ], 0.0),
            MOI.EqualTo(0.0)
        )
    end

    MOI.add_constraint(m, MOI.ScalarAffineFunction([
            MOI.ScalarAffineTerm(1.0, ustart[t])
            for t in 1:L
        ], 0.0),
        MOI.EqualTo(1.0)
    )

    # ==========================================
    #    IV. Add sub-resources to current model
    # ==========================================
    # (None here)
    
    return nothing
end