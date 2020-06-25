"""
    Battery
"""
struct Battery <: Resource
    index::Int           # Index of the resource
    num_timesteps::Int   # Number of time-steps in the battery's operation

    #==================================================
        Battery dynamics
    ==================================================#
    soc_min::Vector{Float64}     # Minimum state of charge for each time step
    soc_max::Vector{Float64}     # Maximum state of charge for each time step
    soc_init::Float64            # Initial state of charge
    selfdischargerate::Float64   # Self-discharging rate
    charge_pwr_min::Float64      # Minimum charging power rate
    charge_pwr_max::Float64      # Maximum charging power rate
    discharge_pwr_min::Float64   # Minimum discharging power rate
    discharge_pwr_max::Float64   # Maximum discharging power rate
    charge_eff::Float64          # Charging efficiency
    discharge_eff::Float64       # Discharging efficiency

    function Battery(;
        index::Integer=0,
        T::Integer=0,
        soc_min::Vector{Float64}=[0.0],
        soc_max::Vector{Float64}=[0.0],
        soc_init::Float64=0.0,
        selfdischargerate::Float64=1.0,
        charge_pwr_min::Float64=0.0,
        charge_pwr_max::Float64=0.0,
        discharge_pwr_min::Float64=0.0,
        discharge_pwr_max::Float64=0.0,
        charge_eff::Float64=0.0,
        discharge_eff::Float64=0.0
    )

        # Dimension checks
        T == size(soc_min, 1) || throw(DimensionMismatch("Invalid soc_min"))
        T == size(soc_max, 1) || throw(DimensionMismatch("Invalid soc_max"))
        0.0 <= charge_eff <= 1.0 || throw(DomainError("charge_eff must be in [0.0, 1.0]"))
        0.0 <= discharge_eff <= 1.0 || throw(DomainError("discharge_eff must be in [0.0, 1.0]"))
        0.0 <= selfdischargerate <= 1.0 || throw(DomainError("selfdischargerate must be in [0.0, 1.0]"))

        return new(
            index, T,
            copy(soc_min), copy(soc_max), soc_init,
            selfdischargerate,
            charge_pwr_min, charge_pwr_max,
            discharge_pwr_min, discharge_pwr_max,
            charge_eff, discharge_eff
        )
    end
end

function add_resource_to_model!(
    m::MOI.ModelLike,
    h::House,
    bat::Battery,
    var2idx, con2idx
)
    T = bat.num_timesteps

    # ==========================================
    #    I. Add local variables
    # ==========================================
    pnet = MOI.add_variables(m, T)  # battery net power
    pchg = MOI.add_variables(m, T)  # battery charging power
    pdis = MOI.add_variables(m, T)  # battery discharging power
    bsoc = MOI.add_variables(m, T)  # State of charge
    uchg = MOI.add_variables(m, T)  # Charge indicator
    udis = MOI.add_variables(m, T)  # Discharge indicator

    # Variable bounds
    for t in 1:T
        MOI.add_constraint(m, MOI.SingleVariable(pchg[t]), MOI.GreaterThan(0.0))
        MOI.add_constraint(m, MOI.SingleVariable(pdis[t]), MOI.GreaterThan(0.0))
        MOI.add_constraint(m, MOI.SingleVariable(bsoc[t]),
            MOI.Interval(bat.soc_min[t], bat.soc_max[t])
        )
        MOI.add_constraint(m, MOI.SingleVariable(uchg[t]), MOI.ZeroOne())
        MOI.add_constraint(m, MOI.SingleVariable(udis[t]), MOI.ZeroOne())
    end

    # Update variable indices
    for t in 1:T
        var2idx[(:bat, bat.index, :pnet, t)] = pnet[t]
        var2idx[(:bat, bat.index, :pchg, t)] = pchg[t]
        var2idx[(:bat, bat.index, :pdis, t)] = pdis[t]
        var2idx[(:bat, bat.index, :bsoc, t)] = bsoc[t]
        var2idx[(:bat, bat.index, :uchg, t)] = uchg[t]
        var2idx[(:bat, bat.index, :udis, t)] = udis[t]
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
        # Net power
        #   `p[t] - pchg[t] + pdis[t] = 0`
        MOI.add_constraint(m, 
            MOI.ScalarAffineFunction([
                MOI.ScalarAffineTerm( 1.0, pnet[t]),
                MOI.ScalarAffineTerm(-1.0, pchg[t]),
                MOI.ScalarAffineTerm( 1.0, pdis[t])
            ], 0.0),
            MOI.EqualTo(0.0)
        )

        # Charging power
        #   `pchg[t] >= uchg[t] * P_chg_min`
        #   `pchg[t] <= uchg[t] * P_chg_max`
        MOI.add_constraint(m, 
            MOI.ScalarAffineFunction([
                MOI.ScalarAffineTerm( 1.0, pchg[t]),
                MOI.ScalarAffineTerm(-bat.charge_pwr_min, uchg[t])
            ], 0.0),
            MOI.GreaterThan(0.0)
        )
        MOI.add_constraint(m, 
            MOI.ScalarAffineFunction([
                MOI.ScalarAffineTerm( 1.0, pchg[t]),
                MOI.ScalarAffineTerm(-bat.charge_pwr_max, uchg[t])
            ], 0.0),
            MOI.LessThan(0.0)
        )

        # Discharging power
        #   `pdis[t] >= udis[t] * P_dis_min`
        #   `pdis[t] <= udis[t] * P_dis_max`
        MOI.add_constraint(m, 
            MOI.ScalarAffineFunction([
                MOI.ScalarAffineTerm( 1.0, pdis[t]),
                MOI.ScalarAffineTerm(-bat.discharge_pwr_min, udis[t])
            ], 0.0),
            MOI.GreaterThan(0.0)
        )
        MOI.add_constraint(m, 
            MOI.ScalarAffineFunction([
                MOI.ScalarAffineTerm( 1.0, pdis[t]),
                MOI.ScalarAffineTerm(-bat.discharge_pwr_max, udis[t])
            ], 0.0),
            MOI.LessThan(0.0)
        )

        # Charge-discharge
        #   `uchg[t] + udis[t] <= 1.0`
        MOI.add_constraint(m, 
            MOI.ScalarAffineFunction([
                MOI.ScalarAffineTerm(1.0, uchg[t]),
                MOI.ScalarAffineTerm(1.0, udis[t])
            ], 0.0),
            MOI.LessThan(1.0)
        )
    end

    # Energy conservation at time t=1
    #   soc[1] - soc_init = bat.charge_eff * pchg[1] - (1.0 / bat.discharge_eff) * pdis[1]
    MOI.add_constraint(m, MOI.ScalarAffineFunction([
            MOI.ScalarAffineTerm(1.0, bsoc[1]),
            MOI.ScalarAffineTerm(-bat.charge_eff, pchg[1]),
            MOI.ScalarAffineTerm(1.0 / bat.discharge_eff, pdis[1]),
        ], 0.0),
        MOI.EqualTo(bat.soc_init * bat.selfdischargerate)
    )
    # Energy conservation at time t>1
    #   soc[t] - soc[t-1] == bat.charge_eff * pchg[t] - (1.0 / bat.discharge_eff) * pdis[t]
    for t in 2:T
        MOI.add_constraint(m, MOI.ScalarAffineFunction([
            MOI.ScalarAffineTerm(1.0, bsoc[t]),
            MOI.ScalarAffineTerm(-bat.selfdischargerate, bsoc[t-1]),
            MOI.ScalarAffineTerm(-bat.charge_eff, pchg[t]),
            MOI.ScalarAffineTerm(1.0 / bat.discharge_eff, pdis[t]),
        ], 0.0),
        MOI.EqualTo(0.0)
    )
    end
    
    # ==========================================
    #    IV. Add sub-resources to current model
    # ==========================================
    # (None here)

    return nothing
end