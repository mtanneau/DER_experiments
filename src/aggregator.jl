import MathProgBase
const MPB = MathProgBase

"""
    Aggregator
"""
mutable struct Aggregator
    num_timesteps::Int   # Number of time steps

    netload_min::Vector{Float64}    # Minimum net load
    netload_max::Vector{Float64}    # Maximum net load

    price_mkt::Vector{Float64}  # Market price for electricity
    resources::Vector{DER.Resource} # List of DERs
end

"""
    buildmodel(...)

Builds a centralized MILP model of the pool.
"""
function buildmodel(a::Aggregator, solver::MPB.AbstractMathProgSolver)

    T = a.num_timesteps

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


    # ==========================================
    #    I. Add local variables
    # ==========================================
    # Net load
    # `pnet[t]`
    for t in 1:T
        var2idx[(:agg, 0, :pnet, t)] = length(var2idx) + 1
    end
    append!(obj, a.price_mkt)
    append!(varlb, a.netload_min)
    append!(varub, a.netload_max)
    append!(vartypes, [:Cont for _ in 1:T])


    # ==========================================
    #    II. Update linking constraints
    # ==========================================
    # (None here)


    # ==========================================
    #    III. Add local constraints
    # ==========================================
    # Linking constraint for total net load
    #   pnet[a, t] - sum_{r} pnet[r, t] = 0.0
    for t in 1:T
        row2idx[(:agg, 0, :link, t)] = length(row2idx) + 1
    end
    append!(rowlb, zeros(T))
    append!(rowub, zeros(T))
    constr_agg = [row2idx[(:agg, 0, :link, t)] for t in 1:T]
    append!(constrI, constr_agg)
    append!(constrJ, [var2idx[(:agg, 0, :pnet, t)] for t in 1:T])
    append!(constrV, ones(T))
    

    # ==========================================
    #    IV. Add sub-resources to current model
    # ==========================================
    for r in a.resources
        DER.addmodel!(r,
            var2idx, obj, varlb, varub, vartypes,
            row2idx, rowlb, rowub, constrI, constrJ, constrV,
            constr_agg
        )
    end


    # ==========================================
    #    V. Instanciate model
    # ==========================================
    numvar = length(var2idx)
    numcon = length(row2idx)
    A = sparse(constrI, constrJ, constrV, numcon, numvar)
    
    m = MPB.LinearQuadraticModel(solver)
    MPB.loadproblem!(m, A, varlb, varub, obj, rowlb, rowub, :Min)
    MPB.setvartype!(m, vartypes)

    return m
end