"""
    ThermalLoad

Thermal load.
"""
struct ThermalLoad <: Resource
    index::Int
    num_timesteps::Int

    temp_min::Vector{Float64}  # Minimum temperature at each time-step
    temp_max::Vector{Float64}  # Maximum temperature at each time-step
    temp_ext::Vector{Float64}  # Outside temperature at each time-step
    temp_init::Float64         # Initial temperature
    pwr_min::Vector{Float64}  # Minimum power when device is on
    pwr_max::Vector{Float64}  # Maximum power when device is on
    C::Float64  # Heat capacity
    η::Float64  # Power efficiency of the heating unit
    μ::Float64  # Thermal conduction coefficient

    binflag::Bool  # Hard/soft on-off flag

    function ThermalLoad(;
        index=0, T=0,
        temp_min=Float64[], temp_max=Float64[], temp_ext=Float64[], temp_init=0.0,
        pwr_min=Float64[], pwr_max=Float64[],
        C=1.0, η=1.0, μ=1.0,
        binflag=true
    )

        # Sanity checks
        C  > 0.0 || throw(DomainError("C  must be positive"))
        T == length(temp_min) || throw(DimensionMismatch(
            "T=$T but temp_min has size $(length(temp_min))"
        ))
        T == length(temp_max) || throw(DimensionMismatch(
            "T=$T but temp_max has size $(length(temp_max))"
        ))
        T == length(temp_ext) || throw(DimensionMismatch(
            "T=$T but temp_ext has size $(length(temp_ext))"
        ))
        T == length(pwr_min) || throw(DimensionMismatch(
            "T=$T but pwr_min has size $(length(pwr_min))"
        ))
        T == length(pwr_max) || throw(DimensionMismatch(
            "T=$T but pwr_max has size $(length(pwr_max))"
        ))


        return new(
            index, T,
            temp_min, temp_max, temp_ext, temp_init,
            pwr_min, pwr_max,
            C, η, μ,
            binflag
        )
    end
end

function add_resource_to_model!(
    m::MOI.ModelLike,
    h::House,
    l::ThermalLoad,
    var2idx, con2idx
)
    T = l.num_timesteps

    # ==========================================
    #    I. Add local variables
    # ==========================================
    pnet = MOI.add_variables(m, T)  # net load
    temp = MOI.add_variables(m, T)  # net load
    uind = MOI.add_variables(m, T)  # net load

    # Variable bounds
    for t in 1:T
        # Net load
        MOI.add_constraint(m, MOI.SingleVariable(pnet[t]),
            MOI.Interval(l.pwr_min[t], l.pwr_max[t])
        )
        # Inside temperature
        MOI.add_constraint(m, MOI.SingleVariable(temp[t]),
            MOI.Interval(l.temp_min[t], l.temp_max[t])
        )
        # On-off indicator
        if l.binflag
            MOI.add_constraint(m, MOI.SingleVariable(uind[t]),
                MOI.ZeroOne()
            )
        else
            MOI.add_constraint(m, MOI.SingleVariable(uind[t]),
                MOI.Interval(0.0, 1.0)
            )
        end
    end
    # Update variable indices
    for t in 1:T
        var2idx[(:thermal, l.index, :pnet, t)] = pnet[t]
        var2idx[(:thermal, l.index, :temp, t)] = temp[t]
        var2idx[(:thermal, l.index, :uind, t)] = uind[t]
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
    for t in 1:T
        # On-off
        # Pmin[t] * u[t] <= p[t], i.e., p[t] - Pmin[t] * u[t] >= 0
        MOI.add_constraint(m,
            MOI.ScalarAffineFunction{Float64}([
                MOI.ScalarAffineTerm(1.0, pnet[t]),
                MOI.ScalarAffineTerm(-l.pwr_min[t], uind[t])
            ], 0.0),
            MOI.GreaterThan(0.0)
        )
        # p[t] <= Pmax[t] * u[t], i.e., p[t] - Pmax[t] * u[t] <= 0
        MOI.add_constraint(m,
            MOI.ScalarAffineFunction{Float64}([
                MOI.ScalarAffineTerm(1.0, pnet[t]),
                MOI.ScalarAffineTerm(-l.pwr_max[t], uind[t])
            ], 0.0),
            MOI.LessThan(0.0)
        )
    end

    # Heat transfer, t = 1
    MOI.add_constraint(m, 
        MOI.ScalarAffineFunction([
            MOI.ScalarAffineTerm((l.η / l.C), pnet[1]),
            MOI.ScalarAffineTerm(-1.0, temp[1])
        ], 0.0),
        MOI.EqualTo(-(l.μ / l.C) * l.temp_ext[1] + (-1.0 + l.μ / l.C) * l.temp_init)
    )

    # Heat transfer, t > 1
    for t in 2:T
        MOI.add_constraint(m, 
        MOI.ScalarAffineFunction([
            MOI.ScalarAffineTerm((l.η / l.C), pnet[t]),
            MOI.ScalarAffineTerm(-1.0, temp[t]),
            MOI.ScalarAffineTerm(1.0 - (l.μ / l.C), temp[t-1])
        ], 0.0),
        MOI.EqualTo(-(l.μ / l.C) * l.temp_ext[t])
    )
    end

    # ==========================================
    #    IV. Add sub-resources to current model
    # ==========================================
    # (None here)

    return nothing
end