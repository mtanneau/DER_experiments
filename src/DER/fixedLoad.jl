"""
    FixedLoad

Fixed load, a.k.a. Uncontrollable load
"""
mutable struct FixedLoad <: Resource
    index::Int           # Index of the resource
    num_timesteps::Int   # Number of time-steps in the battery's operation

    load::Vector{Float64}     # Fixed load at each time step
end

function FixedLoad(;
    index::Integer=0,
    T::Integer=0,
    dt::Float64=1.0,
    load::AbstractVector{T1}=[0.0]
) where{T1<:Real}

    # Dimension checks
    T == size(load, 1) || throw(DimensionMismatch("Invalid soc_min"))
    0.0 < dt || throw(DomainError("dt must be positive"))

    FixedLoad(index, T, load)
end

"""
    LindaOracleMIP(l::FixedLoad, solver::AbstractMathProgSolver)

Construct a MIP oracle that schedules the Battery.
"""
function LindaOracleMIP(l::FixedLoad, solver::MPB.AbstractMathProgSolver)
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
        A_link[t, var2idx[(:fixed, l.index, :pnet, t)]] = 1.0
        A_link[T+t, var2idx[(:fixed, l.index, :pnet, t)]] = 1.0
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
    l::FixedLoad,
    var2idx, obj, varlb, varub, vartypes,
    row2idx, rowlb, rowub, constrI, constrJ, constrV,
    constr_
)
    T = l.num_timesteps
    T_ = length(constr_)
    @assert T_ == 0 || T_ == T

    # Create local variables
    for t in 1:T
        var2idx[(:fixed, l.index, :pnet, t)] = length(var2idx)+1  # net load
    end
    append!(obj, zeros(T))
    append!(varlb, l.load)
    append!(varub, l.load)
    append!(vartypes, [:Cont for t in 1:T])

    # Add local variables to existing linking constraints
    if T_ > 0
        append!(constrI, constr_)
        append!(constrJ, [var2idx[(:fixed, l.index, :pnet, t)] for t in 1:T])
        append!(constrV, -ones(T))
    end

    # Create local constraints
    # (None here)

    # Add sub-resources to current model
    # (none here)

    return nothing
end