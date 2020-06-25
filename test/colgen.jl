import Linda

import GLPK
import MathOptInterface
const MOI = MathOptInterface

import Tulip

# GLPK does not support range constraints, so we bridge
MOI.get(::GLPK.Optimizer, ::MOI.BarrierIterations) = 0
MOI.get(::GLPK.Optimizer, ::MOI.SimplexIterations) = 0


function test_colgen1()
    # Here we optimize a problem with one load in one house
    # Aggregator must enforce more restrictive total load bounds, 
    # and coordinate for price.

    # Problem data
    T = 2
    ptou = [1.0, 1.0]
    pmkt = [1.0, 5.0]
    ptot_min = [0.0, 0.0]
    ptot_max = [1.0, 1.0]

    R = 1  # number of resources

    # Instantiate master problem
    rmp = MOI.Bridges.full_bridge_optimizer(GLPK.Optimizer(), Float64)
    MOI.set(rmp, MOI.Silent(), true)

    fobj = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 0.0)
    MOI.set(rmp, MOI.ObjectiveFunction{typeof(fobj)}(), fobj)
    MOI.set(rmp, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    pagg = MOI.add_variables(rmp, T)  # Aggregator net load

    # Bounds
    for t in 1:T
        MOI.add_constraint(rmp,
            MOI.SingleVariable(pagg[t]),
            MOI.Interval(ptot_min[t], ptot_max[t])
        )
    end

    # Set objective coefficient
    for t in 1:T
        MOI.modify(rmp,
            MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
            MOI.ScalarCoefficientChange(pagg[t], pmkt[t] - ptou[t])
        )
    end

    # Convexity constraints
    con_cvx = MOI.ConstraintIndex[]
    for r in 1:R
        cidx = MOI.add_constraint(rmp,
            MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 0.0),
            MOI.EqualTo(1.0)
        )
        push!(con_cvx, cidx)
    end

    # Linking constraints
    # This is simply -pagg[t] + p1[t] + ... + pR[t] + yp[t] - ym[t] == 0
    con_link = MOI.ConstraintIndex[]
    for t in 1:T
        cidx = MOI.add_constraint(rmp, 
            MOI.ScalarAffineFunction([
                MOI.ScalarAffineTerm(-1.0, pagg[t])
            ], 0.0),
            MOI.EqualTo(0.0)
        )
        push!(con_link, cidx)
    end

    mp = Linda.Master(rmp, con_cvx, con_link, zeros(T), pagg, MOI.VariableIndex[])

    # Create the resource
    l1 = DR.DER.DeferrableLoad(
        index=0, T=T, ncycles=1,
        cycle_begin=[1], cycle_end=[T],
        cycle_energy_min=[1.5], cycle_energy_max=[2.0],
        pwr_min=zeros(T), pwr_max=2.0 * ones(T),
        binflag=false
    )

    l2 = DR.DER.ShiftableLoad(
        T=T,
        cycle_begin_min=1,
        cycle_begin_max=2,
        cycle_length=1,
        cycle_load=[1.0]
    )

    # Create sub-problem
    h1 = DR.DER.House(
        1, T,
        [0.0, 0.0], [2.0, 2.0], ptou,
        [l1, l2]
    )
    sp1 = MOI.Bridges.full_bridge_optimizer(GLPK.Optimizer(), Float64)
    o1 = DR.build_oracle(sp1, h1)

    # Add one artifical column
    # This ensures that the master is feasible
    col0 = Linda.Column(
        100.0,
        h1.index,
        [1, 2], [0.0, 0.0], true
    )
    Linda.add_column!(mp, col0)

    # Solve with column generation
    env = Linda.LindaEnv()
    env.verbose.val = 1
    env.num_cgiter_max.val = 10
    Linda.solve_colgen!(env, mp, [o1])

    # Check solution status
    pagg_ = [
        MOI.get(rmp, MOI.VariablePrimal(), pagg[t])
        for t in 1:T
    ]

    # Solution should be optimal
    @test mp.mp_status == MOI.OPTIMAL

    @test isapprox(pagg_[1], 1.0, atol = 1e-5)
    @test isapprox(pagg_[2], 0.5, atol = 1e-5)

    return nothing
end

function test_colgen2()
    # Here we optimize a model with one load in a first house, and a battery
    # in a second house.
    # The aggregator must coordinate so that the battery is charged at the first
    # period and discharged at the second.

    # Problem data
    T = 2
    # Prices
    # Houses do not see that second period is in fact more expensive
    ptou = [2.0, 1.0]
    pmkt = [1.0, 5.0]
    ptot_min = [0.0, 0.0]
    ptot_max = [10.0, 10.0]

    R = 2  # number of resources

    # Instantiate master problem
    rmp = MOI.Bridges.full_bridge_optimizer(GLPK.Optimizer(), Float64)
    MOI.set(rmp, MOI.Silent(), true)

    fobj = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 0.0)
    MOI.set(rmp, MOI.ObjectiveFunction{typeof(fobj)}(), fobj)
    MOI.set(rmp, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    pagg = MOI.add_variables(rmp, T)  # Aggregator net load
    # We add artificial variables to ensure RMP is feasible
    yp = MOI.add_variables(rmp, T)
    ym = MOI.add_variables(rmp, T)

    # Bounds
    for t in 1:T
        MOI.add_constraint(rmp,
            MOI.SingleVariable(pagg[t]),
            MOI.Interval(ptot_min[t], ptot_max[t])
        )

        MOI.add_constraint(rmp, MOI.SingleVariable(yp[t]), MOI.GreaterThan(0.0))
        MOI.add_constraint(rmp, MOI.SingleVariable(ym[t]), MOI.GreaterThan(0.0))
    end

    # Set objective coefficient
    for t in 1:T
        MOI.modify(rmp,
            MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
            MOI.ScalarCoefficientChange(pagg[t], pmkt[t] - ptou[t])
        )

        MOI.modify(rmp,
            MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
            MOI.ScalarCoefficientChange(yp[t], 1000.0)
        )
        MOI.modify(rmp,
            MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
            MOI.ScalarCoefficientChange(ym[t], 1000.0)
        )
    end

    # Convexity constraints
    con_cvx = MOI.ConstraintIndex[]
    for r in 1:R
        cidx = MOI.add_constraint(rmp,
            MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 0.0),
            MOI.EqualTo(1.0)
        )
        push!(con_cvx, cidx)
    end

    # Linking constraints
    # This is simply -pagg[t] + p1[t] + ... + pR[t] + yp[t] - ym[t] == 0
    con_link = MOI.ConstraintIndex[]
    for t in 1:T
        cidx = MOI.add_constraint(rmp, 
            MOI.ScalarAffineFunction([
                MOI.ScalarAffineTerm(-1.0, pagg[t]),
                MOI.ScalarAffineTerm( 1.0, yp[t]),
                MOI.ScalarAffineTerm(-1.0, ym[t])
            ], 0.0),
            MOI.EqualTo(0.0)
        )
        push!(con_link, cidx)
    end

    mp = Linda.Master(rmp, con_cvx, con_link, zeros(T), pagg, MOI.VariableIndex[])

    # Create the resources
    # House 1: one deferrable load
    l1 = DR.DER.DeferrableLoad(
        index=0, T=T, ncycles=1,
        cycle_begin=[1], cycle_end=[T],
        cycle_energy_min=[3.0], cycle_energy_max=[4.0],
        pwr_min=[0.0, 0.0], pwr_max=[2.0, 2.0],
        binflag=false
    )
    h1 = DR.DER.House(
        1, T,
        [0.0, 0.0], [2.0, 2.0], ptou,
        [l1]
    )
    sp1 = MOI.Bridges.full_bridge_optimizer(GLPK.Optimizer(), Float64)
    o1 = DR.build_oracle(sp1, h1)

    # House 2: one fixed load and one battery
    # House 2 should discharge its battery at T=2
    l2 = DR.DER.FixedLoad(index=0, T=T, load=[1.0, 1.0])
    bat2 = DR.DER.Battery(
        index=0, T=T,
        soc_min=[0.0, 0.0], soc_max=[10.0, 10.0],
        charge_pwr_max=6.0, discharge_pwr_max=6.0,
        charge_eff=1.0,
        discharge_eff=1.0
    )

    # Create sub-problem
    # Second house can inject power into the grid
    h2 = DR.DER.House(
        2, T,
        [-10.0, -10.0], [10.0, 10.0], ptou,
        [l2, bat2]
    )
    sp2 = MOI.Bridges.full_bridge_optimizer(GLPK.Optimizer(), Float64)
    o2 = DR.build_oracle(sp2, h2)

    # Add one initial column per resource
    # This ensures that the master is feasible
    for o in [o1, o2]
        Linda.optimize!(o)
        col, rc = Linda.get_columns(o)[1]
        Linda.add_column!(mp, col)
    end

    # Solve with column generation
    env = Linda.LindaEnv()
    env.verbose.val = 0
    env.num_cgiter_max.val = 10
    Linda.solve_colgen!(env, mp, [o1, o2])

    # Check solution status
    pagg_ = [
        MOI.get(rmp, MOI.VariablePrimal(), pagg[t])
        for t in 1:T
    ]

    # Solution should be optimal
    @test mp.mp_status == MOI.OPTIMAL

    # Check aggregated load
    @test isapprox(pagg_[1], 5.0, atol = 1e-5)
    @test isapprox(pagg_[2], 0.0, atol = 1e-5)

    # TODO: check consumption of each house and battery levels

    return nothing
end

@testset "ColGen" begin
    test_colgen1()
    test_colgen2()
end