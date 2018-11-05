"""
    DeferrableLoad

Deferrable Load
"""
mutable struct DeferrableLoad <: Resource
    index::Int  # Unique index of device
    num_timesteps::Int  # Number of time-steps in horizon
    dt::Float64  # Duration of each time-step

    num_cycles::Int  # Number of cycles
    cycle_begin::Vector{Int}  # Beginning of each cycle
    cycle_end::Vector{Int}    # End of each cycle
    cycle_energy_min::Vector{Float64}  # Minimal energy consumption for each cycle
    cycle_energy_max::Vector{Float64}  # Maximal energy consumption for each cycle

    pwr_min::Vector{Float64}  # Minimum power consumption at each time
    pwr_max::Vector{Float64}  # Maximum power consumption at each time
    binary_flag::Bool  # Indicates whether device is hard or soft on-off
end

function DeferrableLoad(;
    index::Integer=0,
    T::Integer=0,
    dt::Float64=1.0,
    num_cycles::Integer=0,
    cycle_begin::Vector{Int}=Int[],
    cycle_end::Vector{Int}=Int[],
    cycle_energy_min::Vector{Float64}=Float64[],
    cycle_energy_max::Vector{Float64}=Float64[],
    pwr_min::Vector{Float64}=Float64[],
    pwr_max::Vector{Float64}=Float64[],
    binary_flag::Bool=true
)
    # Dimension checks
    dt > 0.0 || throw(DomainError("dt must be positive"))
    T == length(pwr_min) || throw(DimensionMismatch(
        "T=$T but pwr_min has size $(length(pwr_min))"
    ))
    T == length(pwr_max) || throw(DimensionMismatch(
        "T=$T but pwr_max has size $(length(pwr_max))"
    ))
    num_cycles == length(cycle_begin) || throw(DimensionMismatch(
        "N=$num_cycles but cycle_begin has size $(length(cycle_begin))"
    ))
    num_cycles == length(cycle_end) || throw(DimensionMismatch(
        "N=$num_cycles but cycle_end has size $(length(cycle_end))"
    ))
    num_cycles == length(cycle_energy_min) || throw(DimensionMismatch(
        "N=$num_cycles but cycle_energy_max has size $(length(cycle_energy_max))"
    ))
    num_cycles == length(cycle_energy_max) || throw(DimensionMismatch(
        "N=$num_cycles but cycle_energy_max has size $(length(cycle_energy_max))"
    ))

    DeferrableLoad(
        index,
        T, dt,
        num_cycles, cycle_begin, cycle_end,
        cycle_energy_min, cycle_energy_max,
        pwr_min, pwr_max,
        binary_flag
    )
end

"""
    LindaOracleMIP

Instanciate a MIP oracle for said resource.
"""
function LindaOracleMIP(l::DeferrableLoad, solver::MPB.AbstractMathProgSolver)
    T = l.num_timesteps

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
    addmodel!(l,
        var2idx, obj, varlb, varub, vartypes,
        row2idx, rowlb, rowub, constrI, constrJ, constrV,
        Int[]
    )

    numvar = length(var2idx)
    numcon = length(row2idx)

    A_link = spzeros(2*T, numvar)
    for t in 1:T
        A_link[t, var2idx[(:defer, l.index, :pnet, t)]] = 1.0
        A_link[T+t, var2idx[(:defer, l.index, :pnet, t)]] = 1.0
    end

    return LindaOracleMIP(
        l.index,
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
function addmodel!(
    l::DeferrableLoad,
    var2idx, obj, varlb, varub, vartypes,
    row2idx, rowlb, rowub, constrI, constrJ, constrV,
    constr_
)
    T = l.num_timesteps
    T_ = length(constr_)

    (T_ == 0) || (T == T_) || error("T=$T but $T_ linking constraints in input")

    # ==========================================
    #    I. Add local variables
    # ==========================================
    # Net load
    for t in 1:T
        var2idx[(:defer, l.index, :pnet, t)] = length(var2idx) + 1
    end
    append!(obj, zeros(T))
    append!(varlb, l.pwr_min)
    append!(varub, l.pwr_max)
    append!(vartypes, fill(:Cont, T))

    # On-Off indicator
    for t in 1:T
        var2idx[(:defer, l.index, :uind, t)] = length(var2idx) + 1
    end
    append!(obj, zeros(T))
    append!(varlb, zeros(T))
    append!(varub, ones(T))
    if l.binary_flag
        append!(vartypes, fill(:Bin, T))
    else
        append!(vartypes, fill(:Cont, T))
    end

    # ==========================================
    #    II. Update linking constraints
    # ==========================================
    if T_ > 0
        append!(constrI, constr_)
        append!(constrJ, [var2idx[(:defer, l.index, :pnet, t)] for t in 1:T])
        append!(constrV, -ones(T))
    end

    # ==========================================
    #    III. Add local constraints
    # ==========================================
    # Energy consumption for each cycle
    # E_min[k] <= p[t_1] + ... + p[t_k] <= E_max[k]
    for k in 1:l.num_cycles
        # create constraint for cycle k
        i = length(row2idx) + 1
        row2idx[:defer, l.index, :cons, 10000*k + t] = i
        
        L = l.cycle_end[k] - l.cycle_begin[k] + 1  # length of the cycle
        append!(constrI, fill(i, L))
        append!(constrJ, [
            var2idx[(:defer, l.index, :pnet, t)]
            for t in Base.UnitRange(l.cycle_begin[k], l.cycle_end[k])
        ])
        append!(constrV, ones(L))

    end
    append!(rowlb, l.cycle_energy_min)
    append!(rowub, l.cycle_energy_max)

    # On-off
    # Pmin[t] * u[t] <= p[t], i.e.:
    # p[t] - Pmin[t] * u[t] >= 0
    for t in 1:T
        i = length(row2idx) + 1
        row2idx[(:defer, l.index, :pmin, t)] = i

        append!(constrI, [i, i])
        append!(constrJ, [
            var2idx[(:defer, l.index, :pnet, t)],
            var2idx[(:defer, l.index, :uind, t)]
        ])
        append!(constrV, [1.0, -l.pwr_min[t]])
    end
    append!(rowlb, zeros(T))
    append!(rowub, fill(Inf, T))

    # p[t] <= Pmax * u[t], i.e. p[t] - Pmax[t] * u[t] <= 0
    for t in 1:T
        i = length(row2idx) + 1
        row2idx[(:defer, l.index, :pmax, t)] = i
        append!(constrI, [i, i])
        append!(constrJ, [
            var2idx[(:defer, l.index, :pnet, t)],
            var2idx[(:defer, l.index, :uind, t)]
        ])
        append!(constrV, [1.0, -l.pwr_max[t]])
    end
    append!(rowlb, fill(-Inf, T))
    append!(rowub, zeros(T))

    # ==========================================
    #    IV. Add sub-resources to current model
    # ==========================================
    # (None here)
    
    return nothing
end