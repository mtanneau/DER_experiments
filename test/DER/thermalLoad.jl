function test_thermalLoad_1()

    T = 2

    # With these values, the heater uses 0.2 power to keep a ΔT of 1 
    l = DR.DER.ThermalLoad(
        index=1, T=2,
        temp_min=[20.0, 20.0],
        temp_max=[25.0, 25.0],
        temp_ext=[0.0, 0.0],
        temp_init=20.0,
        pwr_min=[0.0, 0.0],
        pwr_max=[10.0, 10.0],
        C=1.0, η=1.0, μ=0.2,
        binflag=false
    )

    # GLPK does not support -in-Interval constraints, so we bridge
    model = MOI.Bridges.full_bridge_optimizer(GLPK.Optimizer(), Float64)

    lmin = [0.0, 0.0]
    lmax = [10.0, 10.0]
    price = [2.0, 1.0]

    # Instantiate initial model
    h = DR.DER.House(
        0, T,
        lmin, lmax, price,
        [l]
    )
    
    model, var2idx, con2idx = DR.DER.build_model!(h, model)

    MOI.optimize!(model)

    @test MOI.get(model, MOI.ObjectiveValue()) ≈ 12.0

    # Temperature should be kept at minimum value, i.e., 20
    T1 = var2idx[(:thermal, l.index, :temp, 1)]
    T2 = var2idx[(:thermal, l.index, :temp, 2)]
    @test MOI.get(model, MOI.VariablePrimal(), T1) ≈ 20.0
    @test MOI.get(model, MOI.VariablePrimal(), T2) ≈ 20.0

    # Check net load
    pnet1 = var2idx[(:thermal, l.index, :pnet, 1)]
    pnet2 = var2idx[(:thermal, l.index, :pnet, 2)]
    @test MOI.get(model, MOI.VariablePrimal(), pnet1) ≈ 4.0
    @test MOI.get(model, MOI.VariablePrimal(), pnet2) ≈ 4.0

    return nothing
    
end

function test_thermalLoad_2()

    T = 2

    # With these values, the heater uses 0.2 power to keep a ΔT of 1 
    l = DR.DER.ThermalLoad(
        index=1, T=2,
        temp_min=[20.0, 20.0],
        temp_max=[25.0, 25.0],
        temp_ext=[0.0, 0.0],
        temp_init=20.0,
        pwr_min=[0.0, 0.0],
        pwr_max=[10.0, 10.0],
        C=1.0, η=1.0, μ=0.2,
        binflag=false
    )

    # GLPK does not support -in-Interval constraints, so we bridge
    model = MOI.Bridges.full_bridge_optimizer(GLPK.Optimizer(), Float64)

    lmin = [0.0, 0.0]
    lmax = [10.0, 10.0]
    price = [1.0, 2.0]  # Now it's worth heating more in the first period

    # Instantiate initial model
    h = DR.DER.House(
        0, T,
        lmin, lmax, price,
        [l]
    )
    
    model, var2idx, con2idx = DR.DER.build_model!(h, model)

    MOI.optimize!(model)

    # T1 should be higher
    # T2 should be at minimum, i.e., 20
    T1 = var2idx[(:thermal, l.index, :temp, 1)]
    T2 = var2idx[(:thermal, l.index, :temp, 2)]
    @test MOI.get(model, MOI.VariablePrimal(), T1) ≈ 25.0
    @test MOI.get(model, MOI.VariablePrimal(), T2) ≈ 20.0

    # Check net load
    pnet1 = var2idx[(:thermal, l.index, :pnet, 1)]
    pnet2 = var2idx[(:thermal, l.index, :pnet, 2)]
    @test MOI.get(model, MOI.VariablePrimal(), pnet1) ≈ 9.0
    @test MOI.get(model, MOI.VariablePrimal(), pnet2) ≈ 0.0

    # Check objective value
    @test MOI.get(model, MOI.ObjectiveValue()) ≈ 9.0
    return nothing
    
end

@testset "Thermal" begin
    test_thermalLoad_1()
    test_thermalLoad_2()
end