"""
    DeferrableLoad

Deferrable Load
"""
struct DeferrableLoad <: Resource
    index::Int  # Unique index of device
    num_timesteps::Int  # Number of time-steps in horizon

    ncycles::Int  # Number of cycles
    cycle_begin::Vector{Int}  # Beginning of each cycle
    cycle_end::Vector{Int}    # End of each cycle
    cycle_energy_min::Vector{Float64}  # Minimal energy consumption for each cycle
    cycle_energy_max::Vector{Float64}  # Maximal energy consumption for each cycle

    pwr_min::Vector{Float64}  # Minimum power consumption at each time
    pwr_max::Vector{Float64}  # Maximum power consumption at each time
    binflag::Bool  # Indicates whether device is hard or soft on-off

    function DeferrableLoad(;
        index::Integer=0,
        T::Integer=0,
        ncycles::Integer=0,
        cycle_begin::Vector{Int}=Int[],
        cycle_end::Vector{Int}=Int[],
        cycle_energy_min::Vector{Float64}=Float64[],
        cycle_energy_max::Vector{Float64}=Float64[],
        pwr_min::Vector{Float64}=Float64[],
        pwr_max::Vector{Float64}=Float64[],
        binflag::Bool=true
    )
        # Dimension checks
        T == length(pwr_min) || throw(DimensionMismatch(
            "T=$T but pwr_min has size $(length(pwr_min))"
        ))
        T == length(pwr_max) || throw(DimensionMismatch(
            "T=$T but pwr_max has size $(length(pwr_max))"
        ))
        ncycles == length(cycle_begin) || throw(DimensionMismatch(
            "N=$ncycles but cycle_begin has size $(length(cycle_begin))"
        ))
        ncycles == length(cycle_end) || throw(DimensionMismatch(
            "N=$ncycles but cycle_end has size $(length(cycle_end))"
        ))
        ncycles == length(cycle_energy_min) || throw(DimensionMismatch(
            "N=$ncycles but cycle_energy_max has size $(length(cycle_energy_max))"
        ))
        ncycles == length(cycle_energy_max) || throw(DimensionMismatch(
            "N=$ncycles but cycle_energy_max has size $(length(cycle_energy_max))"
        ))

        return new(
            index, T,
            ncycles, cycle_begin, cycle_end,
            cycle_energy_min, cycle_energy_max,
            pwr_min, pwr_max,
            binflag
        )
    end
end

"""
    addmodel!(l,
        var2idx, obj, varlb, varub, vartypes,
        row2idx, rowlb, rowub, constrI, constrJ, constrV,
        constr_
    )

Add self to existing model.

# Arguments
- `l`
- `var2idx`
- `obj`
- `varlb`
- `vartypes`
- `row2idx`
- `rowlb`
- `rowub`
- `constrI`
- `constrJ`
- `constrV`
- `constr_`
"""
function add_resource_to_model!(
    m::MOI.ModelLike,
    h::House,
    l::DeferrableLoad,
    var2idx, con2idx
)
    T = l.num_timesteps

    # ==========================================
    #    I. Add local variables
    # ==========================================
    pnet = MOI.add_variables(m, T)  # Net load
    defer = MOI.add_variables(m, T)  # On-Off indicator
    # Update list of variable indices
    for t in 1:T
        var2idx[(:curt, l.index, :pnet, t)] = pnet[t]
        var2idx[(:curt, l.index, :uind, t)] = defer[t]
    end
    
    # Variable bounds
    for t in 1:T
        MOI.add_constraint(m,
            MOI.SingleVariable(pnet[t]),
            MOI.Interval(l.pwr_min[t], l.pwr_max[t])
        )  # min-max load

        if l.binflag
            # Add binary constraint
            MOI.add_constraint(m,
                MOI.SingleVariable(defer[t]),
                MOI.ZeroOne()
            )
        else
            # Just add bounds
            MOI.add_constraint(m,
                MOI.SingleVariable(defer[t]),
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
    # Energy consumption for each cycle
    # E_min[k] <= p[t_1] + ... + p[t_k] <= E_max[k]
    # where t_1, t_k are the first and last period of cycle k
    for k in 1:l.ncycles
        # create constraint for cycle k
        MOI.add_constraint(m,
            MOI.ScalarAffineFunction{Float64}([
                MOI.ScalarAffineTerm(1.0, pnet[t])
                for t in (l.cycle_begin[k]):(l.cycle_end[k])
            ], 0.0),
            MOI.Interval(l.cycle_energy_min[k], l.cycle_energy_max[k])
        )
    end

    # On-off
    for t in 1:T
        # Pmin[t] * u[t] <= p[t], i.e., p[t] - Pmin[t] * u[t] >= 0
        MOI.add_constraint(m,
            MOI.ScalarAffineFunction{Float64}([
                MOI.ScalarAffineTerm(1.0, pnet[t]),
                MOI.ScalarAffineTerm(-l.pwr_min[t], defer[t])
            ], 0.0),
            MOI.GreaterThan(0.0)
        )
        # p[t] <= Pmax * u[t], i.e., p[t] - Pmax[t] * u[t] <= 0
        MOI.add_constraint(m,
            MOI.ScalarAffineFunction{Float64}([
                MOI.ScalarAffineTerm(1.0, pnet[t]),
                MOI.ScalarAffineTerm(-l.pwr_max[t], defer[t])
            ], 0.0),
            MOI.LessThan(0.0)
        )
    end

    # ==========================================
    #    IV. Add sub-resources to current model
    # ==========================================
    # (None here)
    
    return nothing
end