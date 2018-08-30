"""
    Battery
"""
mutable struct Battery <: Resource
    index::Int           # Index of the resource
    num_timesteps::Int   # Number of time-steps in the battery's operation
    dt::Float64          # Length of each time-step

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
end

function Battery(;
    index::Integer=0,
    T::Integer=0,
    dt::Float64=1.0,
    soc_min::AbstractVector{T1}=[0.0],
    soc_max::AbstractVector{T2}=[0.0],
    soc_init::Real=0.0,
    selfdischargerate::Real=1.0,
    charge_pwr_min::Real=0.0,
    charge_pwr_max::Real=0.0,
    discharge_pwr_min::Real=0.0,
    discharge_pwr_max::Real=0.0,
    charge_eff::Real=0.0,
    discharge_eff::Real=0.0
) where{T1<:Real, T2<:Real}

    # Dimension checks
    T == size(soc_min, 1) || throw(DimensionMismatch("Invalid soc_min"))
    T == size(soc_max, 1) || throw(DimensionMismatch("Invalid soc_max"))
    0.0 < dt || throw(DomainError("dt must be positive"))
    0.0 <= charge_eff <= 1.0 || throw(DomainError("charge_eff must be in [0.0, 1.0]"))
    0.0 <= discharge_eff <= 1.0 || throw(DomainError("discharge_eff must be in [0.0, 1.0]"))
    0.0 <= selfdischargerate <= 1.0 || throw(DomainError("selfdischargerate must be in [0.0, 1.0]"))

    Battery(
        index, T, dt,
        copy(soc_min), copy(soc_max), soc_init,
        selfdischargerate,
        charge_pwr_min, charge_pwr_max,
        discharge_pwr_min, discharge_pwr_max,
        charge_eff, discharge_eff
    )
end

"""
    LindaOracleMIP(bat::Battery, solver::AbstractMathProgSolver)

Construct a MIP oracle that schedules the Battery.
"""
function LindaOracleMIP(bat::Battery, solver::MPB.AbstractMathProgSolver)
    T = bat.num_timesteps
    # Construct model
    m = JuMP.Model(solver=solver)

    # add variables
    JuMP.@variable(m, p[1:T])  # net power
    JuMP.@variable(m, pchg[1:T] >= 0.0)  # charging power
    JuMP.@variable(m, pdis[1:T] >= 0.0)  # discharging power
    JuMP.@variable(m, bat.soc_min[t] <= soc[t=1:T] <= bat.soc_max[t])  # state-of-charge
    JuMP.@variable(m, uchg[1:T], Bin)  # charging indicator
    JuMP.@variable(m, udis[1:T], Bin)  # discharging indicator
    
    # set constraints
    for t in 1:T
        # net power balance
        JuMP.@constraint(m, p[t] == pchg[t] - pdis[t])

        # charging power
        JuMP.@constraint(m, uchg[t] * bat.charge_pwr_min <= pchg[t])  # min
        JuMP.@constraint(m, pchg[t] <= uchg[t]*bat.charge_pwr_max)  # max

        # discharging power
        JuMP.@constraint(m, udis[t] * bat.discharge_pwr_min <= pdis[t])  # min
        JuMP.@constraint(m, pdis[t] <= udis[t]*bat.discharge_pwr_max)  # max

        # can't charge and discharge simultaneously
        JuMP.@constraint(m, uchg[t] + udis[t] <= 1)
    end
    
    # energy conservation at time t=0
    JuMP.@constraint(
        m,
        soc[1] == bat.charge_eff * pchg[1] - (1.0 / bat.discharge_eff) * pdis[1]
    )
    for t in 2:T
        # energy conservation at time t>0
        JuMP.@constraint(
            m,
            soc[t] - soc[t-1] == bat.charge_eff * pchg[t] - (1.0 / bat.discharge_eff) * pdis[t]
        )
    end

    JuMP.build(m)
    problem = JuMP.internalmodel(m)


    A_link = spzeros(2*bat.num_timesteps, MPB.numvar(problem))
    for i in 1:bat.num_timesteps
        A_link[i, i] = 1.0
        A_link[bat.num_timesteps+i, i] = 1.0
    end

    return LindaOracleMIP(
        bat.index,
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
    addmodel!(bat::Battery, m::JuMP.Model, constr)

Add battery to existing model.
"""
function addmodel!(bat::Battery, m::JuMP.Model, constr)
    T = bat.num_timesteps

    #========================================
        Local variables
    ========================================#
    # Net power, appears in linking constraints with coefficient `-1.0`
    p = [
        JuMP.@variable(
            m,
            objective=0.0,
            inconstraints=[constr[t]],
            coefficients=[-1.0]
        )
        for t in 1:T
    ]
    # Other local variables
    pchg = JuMP.@variable(m, [1:T], lowerbound = 0.0)  # charging power
    pdis = JuMP.@variable(m, [1:T], lowerbound = 0.0)  # discharging power
    soc = JuMP.@variable(m, [t=1:T], lowerbound=bat.soc_min[t], upperbound=bat.soc_max[t])  # state-of-charge
    uchg = JuMP.@variable(m, [1:T], Bin)  # charging indicator
    udis = JuMP.@variable(m, [1:T], Bin)  # discharging indicator
    

    #========================================
        Local constraints
    ========================================#
    for t in 1:T
        # net power balance
        JuMP.@constraint(m, p[t] == pchg[t] - pdis[t])

        # charging power
        JuMP.@constraint(m, uchg[t] * bat.charge_pwr_min <= pchg[t])  # min
        JuMP.@constraint(m, pchg[t] <= uchg[t]*bat.charge_pwr_max)  # max

        # discharging power
        JuMP.@constraint(m, udis[t] * bat.discharge_pwr_min <= pdis[t])  # min
        JuMP.@constraint(m, pdis[t] <= udis[t]*bat.discharge_pwr_max)  # max

        # can't charge and discharge simultaneously
        JuMP.@constraint(m, uchg[t] + udis[t] <= 1)
    end
    
    # energy conservation at time t=0
    JuMP.@constraint(
        m,
        soc[1] == bat.charge_eff * pchg[1] - (1.0 / bat.discharge_eff) * pdis[1]
    )
    for t in 2:T
        # energy conservation at time t>0
        JuMP.@constraint(
            m,
            soc[t] - soc[t-1] == bat.charge_eff * pchg[t] - (1.0 / bat.discharge_eff) * pdis[t]
        )
    end

    return
end