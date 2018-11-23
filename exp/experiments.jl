import DataFrames
import CSV
using Statistics
using LinearAlgebra
using Random
using BenchmarkTools
using SparseArrays
using Printf

BLAS.set_num_threads(1)

# Demand Response module
include("../src/DemandResponse.jl")
DR = DemandResponse

include("solvers.jl")

"""
    generate_mp(n, T, netloadmin, netloadmax, cols, lpsolver)

Instanciate a Master Problem.

# Arguments
- `num_der`: The number of resources
- `T`: Length of the time-horizon (i.e. number of time-periods)
- `netloadmin`: Minimum total net load
- `netloadmax`: Maximum total net load
- `columns`: Initial set of columns
- `lpsolver`: LP solver for the RMP
"""
function generate_mp(
    num_der, T,
    netloadmin, netloadmax,
    columns,
    lpsolver
)
    # Generate RMP data
    A, obj, rowlb, rowub, collb, colub = generate_RMP_data(
        num_der, T, netloadmin, netloadmax,
        columns, lpsolver
    )

    # Instanciate RMP
    rmp = MPB.LinearQuadraticModel(lpsolver)
    MPB.loadproblem!(rmp, A, collb, colub, obj, rowlb, rowub, :Min)
    
    # Instanciate MP
    mp = Linda.LindaMaster(num_der, 2*T, vcat(netloadmax, netloadmin), 4*T, columns, rmp)
    
    return mp
end

"""
    generate_RMP_data(...)

Generate RMP data
"""
function generate_RMP_data(
    num_der, T, netloadmin, netloadmax,
    columns,
    lpsolver
)

    Alink = vcat(
        spzeros(num_der, 4*T),
        vcat(
            hcat(sparse(1.0I, T, T), sparse(-1.0I, T, T), spzeros(T, 2*T)), # <= u
            hcat(spzeros(T, 2*T), sparse(1.0I, T, T), sparse(-1.0I, T, T)) # >= l
        )
    )
    ncols = length(columns)
    A = hcat(
        Alink,
        vcat(
            hcat([sparsevec([c.idx_subproblem], [1.0], num_der) for c in columns]...),
            hcat([c.col for c in columns]...)
        )
    )
    
    obj = vcat(
        zeros(T),
        fill(1e3, T),
        fill(1e3, T),
        zeros(T),
        [c.cost for c in columns]
    )
    
    rowlb = vcat(ones(num_der), netloadmax, netloadmin)
    rowub = vcat(ones(num_der), netloadmax, netloadmin)
    
    collb = zeros(4*T + ncols)
    colub = fill(Inf, 4*T + ncols)

    return A, obj, rowlb, rowub, collb, colub
end

function generate_RMP_data(
    num_der, T, netloadmin, netloadmax,
    columns,
    lpsolver::Tulip.TulipSolver
)
    ncols = length(columns)
    A = Tulip.TLPLinearAlgebra.UnitBlockAngular(
        Matrix(hcat([c.col for c in columns]...)),
        num_der,
        [c.idx_subproblem for c in columns],
        Matrix(vcat(
            hcat(sparse(1.0I, T, T), sparse(-1.0I, T, T), spzeros(T, 2*T)), # <= u
            hcat(spzeros(T, 2*T), sparse(1.0I, T, T), sparse(-1.0I, T, T)) # >= l
        ))
    )
    
    obj = vcat(
        zeros(T),
        fill(1e3, T),
        fill(1e3, T),
        zeros(T),
        [c.cost for c in columns]
    )
    
    rowlb = vcat(ones(num_der), netloadmax, netloadmin)
    rowub = vcat(ones(num_der), netloadmax, netloadmin)
    
    collb = zeros(4*T + ncols)
    colub = fill(Inf, 4*T + ncols)

    return A, obj, rowlb, rowub, collb, colub
end

"""
    generate_resources()

# Arguments
- `num_der`: The number of DERs
- `T`: The number of time-steps
- `devices_rates`: Ownership rates of each device
- `load_norm`: Normalized reference load profile
- `PV_norm`: Normalized PV output profile
- `price`: Electricity price, in dollar per kW.h.
"""
function generate_resources(
    num_der, T,
    devices_rates::Dict{Symbol, Float64},
    load_norm::Vector,
    PV_norm::Vector,
    temp_ext::Vector,
    price::Vector;
    seed=0
)
    # Sanity checks
    T == length(load_norm) || throw(DimensionMismatch(
        "T=$T but load_norm has length $(length(load_norm))"
    ))
    T == length(PV_norm) || throw(DimensionMismatch(
        "T=$T but load_norm has length $(length(PV_norm))"
    ))
    T == length(price) || throw(DimensionMismatch(
        "T=$T but load_norm has length $(length(price))"
    ))

    Random.seed!(seed)
    resources = Vector{DR.DER.Resource}(undef, num_der)
    
    # Create a pool of `num_der` houses
    # Each household owns appliances according to `devices_rates`
    for i in 1:num_der

        appliances = DR.DER.Resource[]
        γ = 0.5 + rand()  # household's scaling factor

        # ==================
        #    Fixed load     
        # ==================
        # All houses have a fixed load
        if rand() < devices_rates[:fixed]
            l = DR.DER.FixedLoad(
                index=i,
                T=T,
                load= γ .* max.(0.0, (load_norm .+ 0.05 .* randn(T)))
            )
            push!(appliances, l)
        end

        # ==================
        #    Fixed load     
        # ==================
        # Dishwasher
        if rand() < devices_rates[:dishwasher]
            dishwasher = DR.DER.ShiftableLoad(
                index = 10000000 + i,  # TODO: handle index number properly
                T=T,
                cycle_begin_min=1,
                cycle_begin_max=T-1,
                cycle_length=2,
                cycle_load=[1.0, 1.0]
            )
            push!(appliances, dishwasher)
        end

        # ======================
        #    Shiftable load     
        # ======================
        # Clothes washer
        if rand() < devices_rates[:clothes_washer]
            clothes_washer = DR.DER.ShiftableLoad(
                index = 20000000 + i,  # TODO: handle index number properly
                T=T,
                cycle_begin_min=1,
                cycle_begin_max=T-1,
                cycle_length=2,
                cycle_load=[1.0, 1.0]
            )
            push!(appliances, clothes_washer)

            # Clothes dryer
            # Only household with clothes washer may have a clothes dryer
            if rand() < (
                devices_rates[:clothes_dryer] / devices_rates[:clothes_washer]
            )

                clothes_dryer = DR.DER.ShiftableLoad(
                    index = 30000000 + i,  # TODO: handle index number properly
                    T=T,
                    cycle_begin_min=1,
                    cycle_begin_max=T-2,
                    cycle_length=3,
                    cycle_load=[1.0, 1.0, 1.0]
                )
                push!(appliances, clothes_dryer)
            end
        end
        # ===================
        #    Thermal load     
        # ===================
        # Electric heating
        if rand() < devices_rates[:thermal]
            th = DR.DER.ThermalLoad(
                index=i, T=T, dt=1.0,
                temp_min=fill(18.0, T), temp_max=fill(22.0, T),
                temp_init=20.0,
                temp_ext=temp_ext + 0.5 .* randn(T),
                pwr_min=zeros(T), pwr_max=10.0 .* ones(T),
                C=3.0, η=1.0, μ=0.2,
                binary_flag=false
            )
            push!(appliances, th)
        end
        
        # =====================
        #    PV, EV & ES    
        # =====================
        if rand() < devices_rates[:renew]

            # Electric vehicle
            if (T % 24) == 0
                ndays = div(T, 24)

                ev = DR.DER.DeferrableLoad(
                    index=i, T=T, dt=1.0,
                    num_cycles=ndays,
                    cycle_begin=[16 + 24*(d-1) for d in 1:ndays],
                    cycle_end=[24 + 24*(d-1) for d in 1:ndays],
                    cycle_energy_min=fill(10.0, ndays),
                    cycle_energy_max=fill(10.0, ndays),
                    pwr_min=repeat([zeros(15); 1.1*ones(9)], ndays),
                    pwr_max=repeat([zeros(15); 7.7*ones(9)], ndays),
                    binary_flag=true
                )
                push!(appliances, ev)

            else
                # Not an integer number of days
                # Do not include EV
            end

            # PV
            ζ = rand(T)
            pv = DR.DER.CurtailableLoad(
                index=i,
                T=T,
                dt=1.0,
                load= - γ .* PV_norm .* ζ,
                binaryFlag=true
            )
            push!(appliances, pv)

            # Battery
            bat = DR.DER.Battery(
                index=i,
                T=T,
                dt=1.0,
                soc_min=zeros(T),
                soc_max=13.5 .* ones(T),
                soc_init=0.0,
                selfdischargerate=1.0,
                charge_pwr_min=0.0, charge_pwr_max=5.0,
                discharge_pwr_min=0.0, discharge_pwr_max=5.0,
                charge_eff=0.95, discharge_eff=0.95
            )
            push!(appliances, bat)
        end


        # Household
        h = DR.DER.House(
            index=i,
            T=T,
            dt=1.0,
            netload_min = zeros(T),
            netload_max = 10.0 .* ones(T),
            price=price,
            appliances=appliances
        )
        resources[i] = h
    end

    return resources
end

function run_experiment(;
    T0::Int=408+7, T::Int=24, R::Int=1024, ξ::Float64=0.33, seed::Int=42,
    solvers::Vector{Symbol}=Symbol[],
    cg_verbose::Bool=true,
    export_res::Bool=false
)

    # Load data
    ont_load = CSV.read("../dat/load2016.csv")
    ont_price = CSV.read("../dat/price2016.csv")
    ont_prod = CSV.read("../dat/prod2016.csv")
    ont_temp = CSV.read("../dat/temperature2016.csv")

    # Market and Time-Of-Use prices
    p_tou = Vector{Float64}(ont_price.TOU[T0:(T0+T-1)])
    p_mkt = Vector{Float64}(ont_price.HOEP[T0:(T0+T-1)])

    # Normalized load
    load_norm = (ont_load.OntDemand / mean(ont_load.OntDemand))[T0:(T0+T-1)]

    # Normalized PV output
    pv_norm = (ont_prod.SOLAR / mean(ont_prod.SOLAR))[T0:(T0+T-1)]

    # Outside temperature
    temp_ext = ont_temp.Temperature[T0:(T0+T-1)]

    # Min and max total net load
    total_load_min = 0.0 * R * ones(T)
    total_load_max = 6.0 * R * ones(T)

    # Devices' ownership rates
    own_rates = Dict(
        :fixed => 1.0,
        :dishwasher => 0.65,
        :clothes_washer => 0.9,
        :clothes_dryer => 0.75,
        :thermal => 0.6,
        :renew => ξ
    )

    # ===========================
    #   Instanciate resources
    # ===========================
    resources = generate_resources(
        R, T, own_rates,
        load_norm, pv_norm, temp_ext, p_mkt,
        seed=seed
    )

    GRBENV = Gurobi.Env()
    GRBSolverSP = GurobiSolver(
        GRBENV, # Re-use the same Gurobi environment
        OutputFlag=0,
        Threads=1,
        # Tuned parameters
        MIPFocus=1,
        Presolve=-1
    )

    pool = Linda.Oracle.LindaOraclePool(
        [
            Linda.Oracle.LindaOracleMIP(r, GRBSolverSP)
            for r in resources
        ]
    )

    # Generate one initial column for each resource
    initial_cols = Linda.Column[]
    for o in pool.oracles
        Linda.Oracle.query!(o, zeros(2*T), Inf)
        col = Linda.Oracle.get_new_columns(o)[1]
        push!(initial_cols, col)
    end

    # ===========================
    #   Benchmark RMP solvers
    # ===========================

    env = Linda.LindaEnv()
    env[:verbose] = cg_verbose
    env[:num_cgiter_max] = 250
    env[:num_columns_max] = Int(ceil(R/10))

    logs = Dict{Symbol, Any}()

    for s in solvers
        
        # Instanciate MP
        mp_solver = solver(s)
        mp = generate_mp(R, T, total_load_min, total_load_max, initial_cols, mp_solver)

        # Log
        cg_log = Dict{Symbol, Any}(
            :T => T,
            :T0 => T0,
            :R => R,
            :seed => seed,
            :xi => ξ,
            :solver => s
        )

        if cg_verbose
            println("*"^80, "\n")
            println("Problem statistics:")
            println("\tR    = $R")
            println("\tT    = $T")
            println("\tξ    = $ξ")
            println("\tseed = $seed")
            println("\tRMP  = $s")
            println("\tItermax = ", env[:num_cgiter_max])
            println("\tColmax = ", env[:num_columns_max])
            println()
        end

        Linda.solve_colgen!(env, mp, pool, cg_log=cg_log)

        logs[s] = cg_log
        
        # Print output
        if export_res
            df = DataFrames.DataFrame(cg_log)
            CSV.write("res_$(T)_$(R)_$(ξ)_$(seed)_$(s).csv", df)
        end

    end

    return logs

end