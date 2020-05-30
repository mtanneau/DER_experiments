using CSV
using Statistics
using Random

import Gurobi
const GRBENV = Gurobi.Env()

import MathOptInterface
const MOI = MathOptInterface

const SAF = MOI.ScalarAffineFunction
const SAT = MOI.ScalarAffineTerm

import Linda

const DATA_DIR = joinpath(@__DIR__, "../dat")
const RMPDIR = joinpath(@__DIR__, "rmp")

include("../src/DemandResponse.jl")
const DR = DemandResponse

"""
    generate_resources()

# Arguments
- `num_der`: The number of households
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
                index=i, T=T,
                temp_min=fill(18.0, T), temp_max=fill(22.0, T),
                temp_init=20.0,
                temp_ext=temp_ext + 0.5 .* randn(T),
                pwr_min=zeros(T), pwr_max=10.0 .* ones(T),
                C=3.0, η=1.0, μ=0.2,
                binflag=false
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
                    index=i, T=T,
                    ncycles=ndays,
                    cycle_begin=[16 + 24*(d-1) for d in 1:ndays],
                    cycle_end=[24 + 24*(d-1) for d in 1:ndays],
                    cycle_energy_min=fill(10.0, ndays),
                    cycle_energy_max=fill(10.0, ndays),
                    pwr_min=repeat([zeros(15); 1.1*ones(9)], ndays),
                    pwr_max=repeat([zeros(15); 7.7*ones(9)], ndays),
                    binflag=true
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
                load= - γ .* PV_norm .* ζ,
                binflag=true
            )
            push!(appliances, pv)

            # Battery
            bat = DR.DER.Battery(
                index=i,
                T=T,
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
            netload_min = zeros(T),
            netload_max = 10.0 .* ones(T),
            price=price,
            appliances=appliances
        )
        resources[i] = h
    end

    return resources
end

function main(T, R, T0=408+7)

    # Load data
    ont_load = CSV.read(joinpath(DATA_DIR, "load2016.csv"))
    ont_price = CSV.read(joinpath(DATA_DIR, "price2016.csv"))
    ont_prod = CSV.read(joinpath(DATA_DIR, "prod2016.csv"))
    ont_temp = CSV.read(joinpath(DATA_DIR, "temperature2016.csv"))

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
        :renew => 0.33
    )

    # Generate resources
    resources = generate_resources(
        R, T, own_rates,
        load_norm, pv_norm, temp_ext, p_mkt,
        seed=5
    )

    # Instantiate oracles
    oracles = [
    DR.build_oracle(
            MOI.Bridges.full_bridge_optimizer(Gurobi.Optimizer(GRBENV, OutputFlag=0, Threads=1), Float64),
            hh
        )
        for hh in resources
    ]

    # Instantiate the RMP
    rmp = Gurobi.Optimizer(GRBENV);
    MOI.set(rmp, MOI.RawParameter("Threads"), 1)
    MOI.set(rmp, MOI.RawParameter("Method"), 2)
    MOI.set(rmp, MOI.RawParameter("OutputFlag"), 0)
    M = 1e3

    # Linking variables & artificial slacks for linking constraints
    y  = MOI.add_variables(rmp, T)
    zp = MOI.add_variables(rmp, T)
    zm = MOI.add_variables(rmp, T)
    for t in 1:T
        MOI.add_constraint(rmp, MOI.SingleVariable(zp[t]), MOI.GreaterThan(0.0))
        MOI.add_constraint(rmp, MOI.SingleVariable(zm[t]), MOI.GreaterThan(0.0))
        MOI.add_constraint(rmp, MOI.SingleVariable(y[t]), MOI.Interval(total_load_min[t], total_load_max[t]))
    end

    var_link = vcat(y, zp, zm)

    # Linking constraints
    con_link = MOI.add_constraints(rmp,
        [
            SAF([SAT(-1.0, y[t]), SAT(1.0, zp[t]), SAT(-1.0, zm[t])], 0.0)
            for t in 1:T
        ],
        [
            MOI.EqualTo(0.0)
            for t in 1:T
        ]
    );

    # Convexity constraints
    con_cvx = MOI.add_constraints(rmp,
        [SAF(SAT{Float64}[], 0.0) for _ in 1:R],
        [MOI.EqualTo(1.0) for _ in 1:R]
    );

    # Objective
    MOI.set(rmp, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    fobj = SAF{Float64}([
        [
            SAT(M, zp[t]) for t in 1:T
        ]; [
            SAT(M, zm[t]) for t in 1:T
        ]], 0.0
    )
    MOI.set(rmp, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), fobj)

    rhs = zeros(T)
    mp = Linda.Master(rmp, con_cvx, con_link, rhs, var_link, MOI.VariableIndex[])

    # Add initial columns
    for oracle in oracles
        Linda.update!(oracle, false, zeros(T), 0.0)
        Linda.optimize!(oracle)
        
        cols = Linda.get_columns(oracle)
        col, _ = cols[1]
        Linda.add_column!(mp, col)
    end

    env = Linda.LindaEnv()
    env.verbose.val = 1
    env.num_columns_max.val = div(R, 8)
    env.num_cgiter_max.val = 200

    MOI.set(rmp, MOI.Silent(), true)

    cglog = Dict()

    function cg_callback()
        niter::Int = cglog[:n_cg_iter]
        if niter > 0 && niter %10 == 0
            @info "Saving RMP at iter $niter"
            Gurobi.write_model(rmp.inner, joinpath(RMPDIR, "DER_$(T)_$(R)_$(niter).mps.bz2"))
        end
        return nothing
    end

    # Print MP solver
    println("MP solver: ", MOI.get(mp.rmp, MOI.SolverName()), "\n")

    Random.seed!(0)
    Linda.solve_colgen!(env, mp, oracles, cg_log=cglog, callback=cg_callback)
    Gurobi.write_model(rmp.inner, joinpath(RMPDIR, "DER_$(T)_$(R)_$(cglog[:n_cg_iter]).mps.bz2"))
    return
end

if abspath(PROGRAM_FILE) == @__FILE__
    T = parse(Int, ARGS[1])
    R = parse(Int, ARGS[2])
    @info "" T R
    main(T, R)
end