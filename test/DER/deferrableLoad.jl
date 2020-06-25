function test_deferrableLoad()

    T = 2
    load = [2.0, 3.0]

    pmin = zeros(T)
    pmax = ones(T)

    cycle_begin = [1]
    cycle_end = [T]

    emin = [1.2]
    emax = [1.5]

    l = DR.DER.DeferrableLoad(
        index=0, T=T, ncycles=1,
        cycle_begin=cycle_begin, cycle_end=cycle_end,
        cycle_energy_min=emin, cycle_energy_max=emax,
        pwr_min=pmin, pwr_max=pmax,
        binflag=false
    )

    # GLPK does not support -in-Interval constraints, so we bridge
    model = MOI.Bridges.full_bridge_optimizer(GLPK.Optimizer(), Float64)

    lmin = [0.0, 0.0]
    lmax = [10.0, 10.0]
    price = [1.0, -0.01]

    # Instantiate initial model
    h = DR.DER.House(
        0, T,
        lmin, lmax, price,
        [l]
    )

    model, var2idx, con2idx = DR.DER.build_model!(h, model)

    MOI.optimize!(model)

    @test MOI.get(model, MOI.ObjectiveValue()) ≈ 0.2 - 0.01

    # Check net load
    pnet1 = var2idx[(:house, h.index, :pnet, 1)]
    pnet2 = var2idx[(:house, h.index, :pnet, 2)]
    @test MOI.get(model, MOI.VariablePrimal(), pnet1) ≈ 0.2
    @test MOI.get(model, MOI.VariablePrimal(), pnet2) ≈ 1.0

    return nothing
    
end

@testset "Def." begin
    test_deferrableLoad()
end