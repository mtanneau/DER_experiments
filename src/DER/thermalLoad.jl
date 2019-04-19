"""
    ThermalLoad

Thermal load.
"""
struct ThermalLoad <: Resource
    index::Int
    num_timesteps::Int
    dt::Float64

    temp_min::Vector{Float64}  # Minimum temperature at each time-step
    temp_max::Vector{Float64}  # Maximum temperature at each time-step
    temp_ext::Vector{Float64}  # Outside temperature at each time-step
    temp_init::Float64         # Initial temperature
    pwr_min::Vector{Float64}  # Minimum power when device is on
    pwr_max::Vector{Float64}  # Maximum power when device is on
    C::Float64  # Heat capacity
    η::Float64  # Power efficiency of the heating unit
    μ::Float64  # Thermal conduction coefficient

    binary_flag::Bool  # Hard/soft on-off flag
end

function ThermalLoad(;
    index=0, T=0, dt=1.0,
    temp_min=Float64[], temp_max=Float64[], temp_ext=Float64[], temp_init=0.0,
    pwr_min=Float64[], pwr_max=Float64[],
    C=1.0, η=1.0, μ=1.0,
    binary_flag=true
)

    # Sanity checks
    dt > 0.0 || throw(DomainError("dt must be positive"))
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


    ThermalLoad(
        index, T, dt,
        temp_min, temp_max, temp_ext, temp_init,
        pwr_min, pwr_max,
        C, η, μ,
        binary_flag
    )
end

"""
    LindaOracleMIP

Instanciate a MIP oracle for said resource.
"""
function LindaOracleMIP(l::ThermalLoad, solver::MPB.AbstractMathProgSolver)
    T = l.num_timesteps

    var2idx = Dict{Tuple{Symbol, Int, Symbol, Int}, Int}()
    row2idx = Dict{Tuple{Symbol, Int, Symbol, Int}, Int}()
    obj = Float64[]
    varlb = Float64[]
    varub = Float64[]
    vartypes = Symbol[]
    rowlb = Float64[]
    rowub = Float64[]
    constrI = Int[]
    constrJ = Int[]
    constrV = Float64[]

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
        A_link[t, var2idx[(:thermal, l.index, :pnet, t)]] = 1.0
        A_link[T+t, var2idx[(:thermal, l.index, :pnet, t)]] = 1.0
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
    l::ThermalLoad,
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
        var2idx[(:thermal, l.index, :pnet, t)] = length(var2idx) + 1
    end
    append!(obj, zeros(T))
    append!(varlb, l.pwr_min)
    append!(varub, l.pwr_max)
    append!(vartypes, fill(:Cont, T))

    # Inside temperature
    for t in 1:T
        var2idx[(:thermal, l.index, :temp, t)] = length(var2idx) + 1
    end
    append!(obj, zeros(T))
    append!(varlb, l.temp_min)
    append!(varub, l.temp_max)
    append!(vartypes, fill(:Cont, T))

    # On-Off indicator
    for t in 1:T
        var2idx[(:thermal, l.index, :uind, t)] = length(var2idx) + 1
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
        append!(constrJ, [var2idx[(:thermal, l.index, :pnet, t)] for t in 1:T])
        append!(constrV, -ones(T))
    end

    # ==========================================
    #    III. Add local constraints
    # ==========================================
    # Heat transfers at t=1
    # 
    i = length(row2idx) + 1
    row2idx[(:thermal, l.index, :temp, 1)] = i
    push!(rowlb, -(l.μ / l.C) * l.temp_ext[1] + (1.0/l.dt + l.μ / l.C) * l.temp_init)
    push!(rowub, -(l.μ / l.C) * l.temp_ext[1] + (1.0/l.dt + l.μ / l.C) * l.temp_init)
    append!(constrI, [i, i])
    append!(constrJ,[
        var2idx[(:thermal, l.index, :pnet, 1)],
        var2idx[(:thermal, l.index, :temp, 1)]
    ])
    append!(constrV, [
        (l.η / l.C),
        1.0 / l.dt
    ])


    for t in 2:T
        i = length(row2idx) + 1
        row2idx[(:thermal, l.index, :temp, t)] = i
        push!(rowlb, -(l.μ / l.C) * l.temp_ext[t])
        push!(rowub, -(l.μ / l.C) * l.temp_ext[t])
        append!(constrI, [i, i, i])
        append!(constrJ,[
            var2idx[(:thermal, l.index, :pnet, t)],
            var2idx[(:thermal, l.index, :temp, t)],
            var2idx[(:thermal, l.index, :temp, t-1)]
        ])
        append!(constrV, [
            (l.η / l.C),
            1.0 / l.dt,
            - (1.0/l.dt + l.μ / l.C)
        ])
    end

    # On-off
    # Pmin[t] * u[t] <= p[t], i.e.:
    # p[t] - Pmin[t] * u[t] >= 0
    for t in 1:T
        i = length(row2idx) + 1
        row2idx[(:thermal, l.index, :pmin, t)] = i

        append!(constrI, [i, i])
        append!(constrJ, [
            var2idx[(:thermal, l.index, :pnet, t)],
            var2idx[(:thermal, l.index, :uind, t)]
        ])
        append!(constrV, [1.0, -l.pwr_min[t]])
    end
    append!(rowlb, zeros(T))
    append!(rowub, fill(Inf, T))

    # p[t] <= Pmax * u[t], i.e. p[t] - Pmax[t] * u[t] <= 0
    for t in 1:T
        i = length(row2idx) + 1
        row2idx[(:thermal, l.index, :pmax, t)] = i
        append!(constrI, [i, i])
        append!(constrJ, [
            var2idx[(:thermal, l.index, :pnet, t)],
            var2idx[(:thermal, l.index, :uind, t)]
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