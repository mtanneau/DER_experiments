function test_battery()

    T = 2
    bat = DR.DER.Battery(
        index=0, T=T,
        soc_min=zeros(T), soc_max=ones(T),
        charge_pwr_max=1.0,discharge_pwr_max=1.0,
        charge_eff=1.0,
        discharge_eff=1.0
    )

    # GLPK does not support -in-Interval constraints, so we bridge
    model = MOI.Bridges.full_bridge_optimizer(GLPK.Optimizer(), Float64)

    lmin = [-2.0, -2.0]
    lmax = [ 2.0,  2.0]
    price = [1.0, 2.0]

    # Instantiate initial model
    h = DR.DER.House(
        0, T,
        lmin, lmax, price,
        [bat]
    )

    model, var2idx, con2idx = DR.DER.build_model!(h, model)

    MOI.optimize!(model)

    @test MOI.get(model, MOI.ObjectiveValue()) ≈ -1.0

    # Check net load
    # Battery should charge then discharge
    pnet1 = var2idx[(:house, h.index, :pnet, 1)]
    pnet2 = var2idx[(:house, h.index, :pnet, 2)]
    @test MOI.get(model, MOI.VariablePrimal(), pnet1) ≈ 1.0
    @test MOI.get(model, MOI.VariablePrimal(), pnet2) ≈ -1.0

    return nothing
end

@testset "Battery" begin
    test_battery()
end