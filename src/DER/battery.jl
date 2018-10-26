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

    var2idx = Dict{Tuple{Symbol, Int, Symbol, Int}, Int}()
    row2idx = Dict{Tuple{Symbol, Int, Symbol, Int}, Int}()
    obj = Vector{Float64}(undef, 0)
    varlb = Vector{Float64}(undef, 0)
    varub = Vector{Float64}(undef, 0)
    vartypes = Vector{Symbol}(undef, 0)
    rowlb = Vector{Float64}(undef, 0)
    rowub = Vector{Float64}(undef, 0)
    constrI = Vector{Int}(undef, 0)
    constrJ = Vector{Int}(undef, 0)
    constrV = Vector{Float64}(undef, 0)

    # Update model
    addmodel!(bat,
        var2idx, obj, varlb, varub, vartypes,
        row2idx, rowlb, rowub, constrI, constrJ, constrV,
        Vector{Int}(undef, 0)
    )
    
    numvar = length(var2idx)
    numcon = length(row2idx)

    A_link = spzeros(2*T, numvar)
    for t in 1:T
        A_link[t, var2idx[(:bat, bat.index, :pnet, t)]] = 1.0
        A_link[T+t, var2idx[(:bat, bat.index, :pnet, t)]] = 1.0
    end

    return LindaOracleMIP(
        bat.index,
        obj,
        A_link,
        sparse(constrI, constrJ, constrV, numcon, numvar),
        rowlb,
        rowub,
        vartypes,
        varlb,
        varub,
        solver
    )

end

function addmodel!(
    bat::Battery,
    var2idx, obj, varlb, varub, vartypes,
    row2idx, rowlb, rowub, constrI, constrJ, constrV,
    constr_
)
    T = bat.num_timesteps
    T_ = length(constr_)
    @assert T_ == 0 || T_ == T

    # ==========================================
    #    I. Add local variables
    # ==========================================
    numvar = 6*T
    append!(varlb, zeros(numvar))
    append!(varub, zeros(numvar))
    append!(obj, zeros(numvar))
    append!(vartypes, Vector{Symbol}(undef, numvar))

    for t in 1:T
        # net load
        var2idx[(:bat, bat.index, :pnet, t)] = length(var2idx)+1
        varlb[var2idx[(:bat, bat.index, :pnet, t)]] = -Inf
        varub[var2idx[(:bat, bat.index, :pnet, t)]] = Inf
        vartypes[var2idx[(:bat, bat.index, :pnet, t)]] = :Cont

        # Charging power
        var2idx[(:bat, bat.index, :pchg, t)] = length(var2idx)+1
        varlb[var2idx[(:bat, bat.index, :pchg, t)]] = 0.0
        varub[var2idx[(:bat, bat.index, :pchg, t)]] = Inf
        vartypes[var2idx[(:bat, bat.index, :pchg, t)]] = :Cont
        
        # Discharging power
        var2idx[(:bat, bat.index, :pdis, t)] = length(var2idx)+1
        varlb[var2idx[(:bat, bat.index, :pdis, t)]] = 0.0
        varub[var2idx[(:bat, bat.index, :pdis, t)]] = Inf
        vartypes[var2idx[(:bat, bat.index, :pdis, t)]] = :Cont

        # State of charge
        var2idx[(:bat, bat.index, :soc, t)] = length(var2idx)+1
        varlb[var2idx[(:bat, bat.index, :soc, t)]] = bat.soc_min[t]
        varub[var2idx[(:bat, bat.index, :soc, t)]] = bat.soc_max[t]
        vartypes[var2idx[(:bat, bat.index, :soc, t)]] = :Cont

        # Charge indicator
        var2idx[(:bat, bat.index, :uchg, t)] = length(var2idx)+1
        varlb[var2idx[(:bat, bat.index, :uchg, t)]] = 0.0
        varub[var2idx[(:bat, bat.index, :uchg, t)]] = 1.0
        vartypes[var2idx[(:bat, bat.index, :uchg, t)]] = :Bin

        # Discharge indicator
        var2idx[(:bat, bat.index, :udis, t)] = length(var2idx)+1
        varlb[var2idx[(:bat, bat.index, :udis, t)]] = 0.0
        varub[var2idx[(:bat, bat.index, :udis, t)]] = 1.0
        vartypes[var2idx[(:bat, bat.index, :udis, t)]] = :Bin
    end

    # ==========================================
    #    II. Update linking constraints
    # ==========================================
    if T_ > 0
        append!(constrI, constr_)
        append!(constrJ, [var2idx[(:bat, bat.index, :pnet, t)] for t in 1:T])
        append!(constrV, -ones(T))
    end

    # ==========================================
    #    III. Add local constraints
    # ==========================================
    numcon = 7*T
    append!(rowlb, zeros(numcon))
    append!(rowub, zeros(numcon))

    for t in 1:T
        # Net power
        #   `p[t] - pchg[t] + pdis[t] = 0`
        row2idx[(:bat, bat.index, :netload, t)] = length(row2idx)+1
        rowlb[row2idx[(:bat, bat.index, :netload, t)]] = 0.0
        rowub[row2idx[(:bat, bat.index, :netload, t)]] = 0.0
        i = row2idx[(:bat, bat.index, :netload, t)]
        append!(constrI, [i, i, i])
        append!(constrJ,
            [
                var2idx[(:bat, bat.index, :pnet, t)],
                var2idx[(:bat, bat.index, :pchg, t)],
                var2idx[(:bat, bat.index, :pdis, t)]
            ]
        )
        append!(constrV, [1.0, -1.0, 1.0])

        # Charging power
        #   `pchg[t] >= uchg[t] * P_chg_min`
        #   `pchg[t] <= uchg[t] * P_chg_max`
        row2idx[(:bat, bat.index, :pchg_min, t)] = length(row2idx)+1
        rowlb[row2idx[(:bat, bat.index, :pchg_min, t)]] = 0.0
        rowub[row2idx[(:bat, bat.index, :pchg_min, t)]] = Inf
        i = row2idx[(:bat, bat.index, :pchg_min, t)]
        append!(constrI, [i, i])
        append!(constrJ,
            [
                var2idx[(:bat, bat.index, :pchg, t)],
                var2idx[(:bat, bat.index, :uchg, t)]
            ]
        )
        append!(constrV, [1.0, -bat.charge_pwr_min])

        row2idx[(:bat, bat.index, :pchg_max, t)] = length(row2idx)+1
        rowlb[row2idx[(:bat, bat.index, :pchg_max, t)]] = -Inf
        rowub[row2idx[(:bat, bat.index, :pchg_max, t)]] = 0.0
        i = row2idx[(:bat, bat.index, :pchg_max, t)]
        append!(constrI, [i, i])
        append!(constrJ,
            [
                var2idx[(:bat, bat.index, :pchg, t)],
                var2idx[(:bat, bat.index, :uchg, t)]
            ]
        )
        append!(constrV, [1.0, -bat.charge_pwr_max])
        
        # Discharging power
        #   `pdis[t] >= udis[t] * P_dis_min`
        #   `pdis[t] <= udis[t] * P_dis_max`
        row2idx[(:bat, bat.index, :pdis_min, t)] = length(row2idx)+1
        rowlb[row2idx[(:bat, bat.index, :pdis_min, t)]] = 0.0
        rowub[row2idx[(:bat, bat.index, :pdis_min, t)]] = Inf
        i = row2idx[(:bat, bat.index, :pdis_min, t)]
        append!(constrI, [i, i])
        append!(constrJ,
            [
                var2idx[(:bat, bat.index, :pdis, t)],
                var2idx[(:bat, bat.index, :udis, t)]
            ]
        )
        append!(constrV, [1.0, -bat.discharge_pwr_min])

        row2idx[(:bat, bat.index, :pdis_max, t)] = length(row2idx)+1
        rowlb[row2idx[(:bat, bat.index, :pdis_max, t)]] = -Inf
        rowub[row2idx[(:bat, bat.index, :pdis_max, t)]] = 0.0
        i = row2idx[(:bat, bat.index, :pdis_max, t)]
        append!(constrI, [i, i])
        append!(constrJ,
            [
                var2idx[(:bat, bat.index, :pdis, t)],
                var2idx[(:bat, bat.index, :udis, t)]
            ]
        )
        append!(constrV, [1.0, -bat.discharge_pwr_max])


        # Charge-discharge
        #   `uchg[t] + udis[t] <= 1.0`
        row2idx[(:bat, bat.index, :chgdis, t)] = length(row2idx)+1
        rowlb[row2idx[(:bat, bat.index, :chgdis, t)]] = -Inf
        rowub[row2idx[(:bat, bat.index, :chgdis, t)]] = 1.0
        i = row2idx[(:bat, bat.index, :chgdis, t)]
        append!(constrI, [i, i])
        append!(constrJ,
            [
                var2idx[(:bat, bat.index, :uchg, t)],
                var2idx[(:bat, bat.index, :udis, t)]
            ]
        )
        append!(constrV, [1.0, 1.0])
        
    end

    # Energy conservation at time t=1
    #   soc[1] - soc_init = bat.charge_eff * pchg[1] - (1.0 / bat.discharge_eff) * pdis[1]
    row2idx[(:bat, bat.index, :soc, 1)] = length(row2idx)+1
    rowlb[row2idx[(:bat, bat.index, :soc, 1)]] = bat.soc_init
    rowub[row2idx[(:bat, bat.index, :soc, 1)]] = bat.soc_init
    i = row2idx[(:bat, bat.index, :soc, 1)]
    append!(constrI, [i, i, i])
    append!(constrJ,
        [
            var2idx[(:bat, bat.index, :soc, 1)],
            var2idx[(:bat, bat.index, :pchg, 1)],
            var2idx[(:bat, bat.index, :pdis, 1)]
        ]
    )
    append!(constrV,
        [
            1.0,
            -bat.charge_eff,
            1.0 / bat.discharge_eff
        ]
    )

    # Energy conservation at time t>1
    #   soc[t] - soc[t-1] == bat.charge_eff * pchg[t] - (1.0 / bat.discharge_eff) * pdis[t]
    for t in 2:T
        row2idx[(:bat, bat.index, :soc, t)] = length(row2idx)+1
        rowlb[row2idx[(:bat, bat.index, :soc, t)]] = 0.0
        rowub[row2idx[(:bat, bat.index, :soc, t)]] = 0.0
        i = row2idx[(:bat, bat.index, :soc, t)]
        append!(constrI, [i, i, i, i])
        append!(constrJ,
            [
                var2idx[(:bat, bat.index, :soc, t)],
                var2idx[(:bat, bat.index, :soc, t-1)],
                var2idx[(:bat, bat.index, :pchg, t)],
                var2idx[(:bat, bat.index, :pdis, t)]
            ]
        )
        append!(constrV,
            [
                1.0,
                -1.0,
                -bat.charge_eff,
                1.0 / bat.discharge_eff
            ]
        )
    end

    # ==========================================
    #    IV. Add sub-resources to current model
    # ==========================================
    # (None here)

    return nothing
end