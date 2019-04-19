"""
    ShiftableLoad

ShiftableLoad Load
"""
struct ShiftableLoad <: Resource
    index::Int  # Unique index of device
    num_timesteps::Int  # Number of time-steps in horizon
    dt::Float64  # Duration of each time-step

    cycle_begin_min::Int  # Beginning of each cycle
    cycle_begin_max::Int    # End of each cycle
    cycle_length::Int  # Length of a cycle
    cycle_load::Vector{Float64}  # Load level at each time of the cycle
end

function ShiftableLoad(;
    index::Integer=0,
    T::Integer=0,
    dt::Float64=1.0,
    cycle_begin_min::Int=1,
    cycle_begin_max::Int=0,
    cycle_length::Int=0,
    cycle_load::Vector{Float64}=Float64[]
)
    # Dimension checks
    dt > 0.0 || throw(DomainError("dt must be positive"))
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

    ShiftableLoad(
        index,
        T, dt,
        cycle_begin_min, cycle_begin_max,
        cycle_length, cycle_load
    )
end

"""
    LindaOracleMIP

Instanciate a MIP oracle for said resource.
"""
function LindaOracleMIP(l::ShiftableLoad, solver::MPB.AbstractMathProgSolver)
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
        A_link[t, var2idx[(:shift, l.index, :pnet, t)]] = 1.0
        A_link[T+t, var2idx[(:shift, l.index, :pnet, t)]] = 1.0
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
    l::ShiftableLoad,
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
        var2idx[(:shift, l.index, :pnet, t)] = length(var2idx) + 1
    end
    append!(obj, zeros(T))
    append!(varlb, fill(-Inf, T))
    append!(varub, fill(Inf, T))
    append!(vartypes, fill(:Cont, T))

    # Cycle start indicator
    for t in l.cycle_begin_min:l.cycle_begin_max
        var2idx[(:shift, l.index, :start, t)] = length(var2idx) + 1
    end
    append!(obj, zeros(l.cycle_begin_max - l.cycle_begin_min + 1))
    append!(varlb, zeros(l.cycle_begin_max - l.cycle_begin_min + 1))
    append!(varub, ones(l.cycle_begin_max - l.cycle_begin_min + 1))
    append!(vartypes, fill(:Bin, l.cycle_begin_max - l.cycle_begin_min + 1))

    # ==========================================
    #    II. Update linking constraints
    # ==========================================
    if T_ > 0
        append!(constrI, constr_)
        append!(constrJ, [var2idx[(:shift, l.index, :pnet, t)] for t in 1:T])
        append!(constrV, -ones(T))
    end

    # ==========================================
    #    III. Add local constraints
    # ==========================================
    # Cycle must start-up exactly once
    # v[t1]+ ... + v[tL] = 1
    i = length(row2idx) + 1
    row2idx[(:shift, l.index, :start, 1)] = i
    push!(rowlb, 1.0)
    push!(rowub, 1.0)
    append!(constrI, [i   for _ in l.cycle_begin_min:l.cycle_begin_max])
    append!(constrJ, [
        var2idx[(:shift, l.index, :start, t)]
        for t in l.cycle_begin_min:l.cycle_begin_max
    ])
    append!(constrV, [1.0 for _ in l.cycle_begin_min:l.cycle_begin_max])

    # Net load
    for t in 1:T
        i = length(row2idx) + 1
        row2idx[(:shift, l.index, :pnet, t)] = i

        append!(constrI, [
            [i];
            [   
                i for k in 1:l.cycle_length
                if l.cycle_begin_min <= t-k+1 <= l.cycle_begin_max
            ]
        ])
        append!(constrJ, [
            [var2idx[(:shift, l.index, :pnet, t)]];
            [
                var2idx[(:shift, l.index, :start, t-k+1)]
                for k in 1:l.cycle_length
                if l.cycle_begin_min <= t-k+1 <= l.cycle_begin_max
            ]
        ])
        append!(constrV, [
            [-1.0];
            [
                l.cycle_load[k]
                for k in 1:l.cycle_length
                if l.cycle_begin_min <= t-k+1 <= l.cycle_begin_max
            ]
        ])

    end
    append!(rowlb, zeros(T))
    append!(rowub, zeros(T))

    # ==========================================
    #    IV. Add sub-resources to current model
    # ==========================================
    # (None here)
    
    return nothing
end