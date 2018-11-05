"""
    CurtailableLoad

Curtailable load
"""
struct CurtailableLoad <: Resource
    index::Int           # Index of the resource
    num_timesteps::Int   # Number of time-steps in the battery's operation

    load::Vector{Float64}     # Load at each time-step, if device is on
    binaryFlag::Bool  # Whether on-off or continuous curtailment
end

function CurtailableLoad(;
    index::Integer=0,
    T::Integer=0,
    dt::Float64=1.0,
    load::Vector{Float64}=Float64[],
    binaryFlag::Bool=true
)

    # Dimension checks
    T == size(load, 1) || throw(DimensionMismatch("Invalid load dimension"))
    0.0 < dt || throw(DomainError("dt must be positive"))

    CurtailableLoad(index, T, load, binaryFlag)
end

"""
    LindaOracleMIP(l::CurtailableLoad, solver::AbstractMathProgSolver)

Construct a MIP oracle that schedules the Battery.
"""
function LindaOracleMIP(l::CurtailableLoad, solver::MPB.AbstractMathProgSolver)
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
        Vector{Int}(undef, 0)
    )

    numvar = length(var2idx)
    numcon = length(row2idx)

    A_link = spzeros(2*T, numvar)
    for t in 1:T
        A_link[t, var2idx[(:curt, l.index, :pnet, t)]] = 1.0
        A_link[T+t, var2idx[(:curt, l.index, :pnet, t)]] = 1.0
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


function addmodel!(
    l::CurtailableLoad,
    var2idx, obj, varlb, varub, vartypes,
    row2idx, rowlb, rowub, constrI, constrJ, constrV,
    constr_
)
    T = l.num_timesteps
    T_ = length(constr_)
    (T_ == 0) || (T_ == T) || error("T=$T but $T_ linking constraints in input")

    
    # ==========================================
    #    I. Add local variables
    # ==========================================
    # Net load
    for t in 1:T
        var2idx[(:curt, l.index, :pnet, t)] = length(var2idx)+1
    end
    append!(obj, zeros(T))
    append!(varlb, fill(-Inf, T))
    append!(varub, fill(Inf, T))
    append!(vartypes, fill(:Cont, T))

    # Curtailment indicator
    for t in 1:T
        var2idx[(:curt, l.index, :uind, t)] = length(var2idx)+1
    end
    append!(obj, zeros(T))
    append!(varlb, zeros(T))
    append!(varub, ones(T))
    if l.binaryFlag
        append!(vartypes, fill(:Bin, T))
    else
        append!(vartypes, fill(:Cont, T))
    end


    # ==========================================
    #    II. Update linking constraints
    # ==========================================
    if T_ > 0
        append!(constrI, constr_)
        append!(constrJ, [var2idx[(:curt, l.index, :pnet, t)] for t in 1:T])
        append!(constrV, -ones(T))
    end

    # ==========================================
    #    III. Add local constraints
    # ==========================================
    # Curtailement
    #   p[t] - u[t] * P[t] = 0
    for t in 1:T
        row2idx[(:curt, l.index, :curtail, t)] = length(row2idx) + 1
        i = row2idx[(:curt, l.index, :curtail, t)]
        append!(constrI, [i, i])
        append!(constrJ, 
            [
                var2idx[(:curt, l.index, :pnet, t)],
                var2idx[(:curt, l.index, :uind, t)]
            ]
        )
        append!(constrV, [1.0, -l.load[t]])
    end
    append!(rowlb, zeros(T))
    append!(rowub, zeros(T))
    

    # ==========================================
    #    IV. Add sub-resources to current model
    # ==========================================
    # (None here)

    return nothing
end