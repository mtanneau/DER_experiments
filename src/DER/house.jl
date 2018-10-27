"""
    House

Data stucture that represents a house and its appliances.

# Attributes
- `index`: Index of that house
- `num_timesteps`: Number of time-steps in the time horizon
- `dt`: Duration of one time-step, in hours
- `netload_min`: House's minimum net load, for each time-step
- `netload_max`:
- `price`: Electricity price, in dollar per kWh
- `appliances`: List of appliances in that house
"""
mutable struct House <: Resource
    index::Int

    num_timesteps::Int
    dt::Float64

    netload_min::Vector{Float64}
    netload_max::Vector{Float64}
    price::Vector{Float64}

    appliances::Vector{Resource}
end

function House(;
    index::Integer=0,
    T::Integer=0,
    dt::Float64=1.0,
    netload_min::Vector{Float64}=Vector{Float64}(),
    netload_max::Vector{Float64}=Vector{Float64}(),
    price::Vector{Float64}=Vector{Float64}(),
    appliances::Vector{Resource}=Vector{Resource}()
)

    # Dimension checks
    T == size(netload_min, 1) || throw(DimensionMismatch("Invalid netload_min"))
    T == size(netload_max, 1) || throw(DimensionMismatch("Invalid netload_max"))
    T == size(price, 1) || throw(DimensionMismatch("Invalid price"))
    0.0 < dt || throw(DomainError("dt must be positive"))

    House(
        index, T, dt,
        copy(netload_min), copy(netload_max), copy(price),
        appliances
    )
end

"""
    LindaOracleMIP(h, solver)

Construct a MIP oracle that schedules the House.
"""
function LindaOracleMIP(h::House, solver::MPB.AbstractMathProgSolver)

    T = h.num_timesteps

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
    addmodel!(h,
        var2idx, obj, varlb, varub, vartypes,
        row2idx, rowlb, rowub, constrI, constrJ, constrV,
        Vector{Int}(undef, 0)
    )
    
    numvar = length(var2idx)
    numcon = length(row2idx)

    A_link = spzeros(2*T, numvar)
    for t in 1:T
        A_link[t, var2idx[(:house, h.index, :pnet, t)]] = 1.0
        A_link[T+t, var2idx[(:house, h.index, :pnet, t)]] = 1.0
    end

    return LindaOracleMIP(
        h.index,
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
    h::House,
    var2idx, obj, varlb, varub, vartypes,
    row2idx, rowlb, rowub, constrI, constrJ, constrV,
    constr_
)
    T = h.num_timesteps
    T_ = length(constr_)
    @assert T_ == 0 || T_ == T

    # ==========================================
    #    I. Add local variables
    # ==========================================
    for t in 1:T
        var2idx[(:house, h.index, :pnet, t)] = length(var2idx)+1  # net load
    end
    append!(obj, h.price)
    append!(varlb, h.netload_min)
    append!(varub, h.netload_max)
    append!(vartypes, [:Cont for t in 1:T])

    # ==========================================
    #    II. Update linking constraints
    # ==========================================
    if T_ > 0
        append!(constrI, constr_)
        append!(constrJ, [var2idx[(:house, h.index, :pnet, t)] for t in 1:T])
        append!(constrV, -ones(T))
    end

    # ==========================================
    #    III. Add local constraints
    # ==========================================

    # Linking constraint for house's net load
    #   pnet[h, t] - sum_{r} pnet[r, t] = 0.0
    for t in 1:T
        row2idx[(:house, h.index, :link, t)] = length(row2idx) + 1
    end
    append!(rowlb, zeros(T))
    append!(rowub, zeros(T))

    constr_house = [row2idx[(:house, h.index, :link, t)] for t in 1:T]
    append!(constrI, constr_house)
    append!(constrJ, [var2idx[(:house, h.index, :pnet, t)] for t in 1:T])
    append!(constrV, ones(T))
    

    # ==========================================
    #    IV. Add sub-resources to current model
    # ==========================================
    for r in h.appliances
        addmodel!(r,
            var2idx, obj, varlb, varub, vartypes,
            row2idx, rowlb, rowub, constrI, constrJ, constrV,
            constr_house
        )
    end

    return nothing
end