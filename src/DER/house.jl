"""
    House

Data stucture that represents a house and its appliances.

# Attributes
- `index`: Index of that house
- `num_timesteps`: Number of time-steps in the time horizon
- `netload_min`: House's minimum net load, for each time-step
- `netload_max`:
- `price`: Electricity price, in dollar per kWh
- `appliances`: List of appliances in that house
"""
struct House <: Resource
    index::Int
    num_timesteps::Int

    netload_min::Vector{Float64}
    netload_max::Vector{Float64}
    price::Vector{Float64}

    appliances::Vector{Resource}
end

function House(;
    index::Integer=0,
    T::Integer=0,
    netload_min::Vector{Float64}=Vector{Float64}(),
    netload_max::Vector{Float64}=Vector{Float64}(),
    price::Vector{Float64}=Vector{Float64}(),
    appliances::Vector{Resource}=Vector{Resource}()
)

    # Dimension checks
    T == size(netload_min, 1) || throw(DimensionMismatch("Invalid netload_min"))
    T == size(netload_max, 1) || throw(DimensionMismatch("Invalid netload_max"))
    T == size(price, 1) || throw(DimensionMismatch("Invalid price"))

    House(
        index, T,
        copy(netload_min), copy(netload_max), copy(price),
        appliances
    )
end

"""
    LindaOracleMIP(h, solver)

Construct a MIP oracle that schedules the House.
"""
function build_model!(h::House, model::MOI.ModelLike)
    # Emtpy model
    # TODO: use an UninstanitatedOptimizer instead
    MOI.empty!(model)

    # Minimize
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    # Initial objective function is empty
    fobj = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 0.0)
    MOI.set(model, MOI.ObjectiveFunction{typeof(fobj)}(), fobj)
    # TODO: set parameters

    T = h.num_timesteps

    # Initialize
    var2idx = Dict{Tuple{Symbol, Int, Symbol, Int}, MOI.VariableIndex}()
    con2idx = Dict{Tuple{Symbol, Int, Symbol, Int}, MOI.ConstraintIndex}()

    # ==========================================
    #    I. Add local variables
    # ==========================================
    pnet = MOI.add_variables(model, T)  # House's net load
    # Variable bounds
    for t in 1:T
        MOI.add_constraint(model,
            MOI.SingleVariable(pnet[t]),
            MOI.Interval(h.netload_min[t], h.netload_max[t])
        )
    end

    # Set objective coefficients (?)
    # append!(obj, h.price)
    for t in 1:T
        MOI.modify(model,
            MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
            MOI.ScalarCoefficientChange(pnet[t], h.price[t])
        )
    end

    # Update list of variable indices
    for t in 1:T
        var2idx[(:house, h.index, :pnet, t)] = pnet[t]  # net load
    end

    # ==========================================
    #    II. Add local constraints
    # ==========================================
    # Linking constraint for house's net load
    #   pnet[h, t] - sum_{r} pnet[r, t] = 0.0
    for t in 1:T
        cidx = MOI.add_constraint(model,
            MOI.ScalarAffineFunction{Float64}(
                [MOI.ScalarAffineTerm(1.0, pnet[t])],
                0.0
            ),
            MOI.EqualTo(0.0)
        )

        # Record constraint index
        con2idx[(:house, h.index, :link, t)] = cidx
    end

    # ==========================================
    #    III. Add sub-resources to current model
    # ==========================================
    for r in h.appliances
        add_resource_to_model!(model, h, r, var2idx, con2idx)
    end

    

    return model, var2idx, con2idx
end