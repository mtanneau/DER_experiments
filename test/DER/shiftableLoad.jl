function test_shiftableLoad()

    T = 4

    l = DR.DER.ShiftableLoad(
        T=T,
        cycle_begin_min=2,
        cycle_begin_max=3,
        cycle_length=2,
        cycle_load=[1.0, 2.0]
    )

    # GLPK does not support -in-Interval constraints, so we bridge
    model = MOI.Bridges.full_bridge_optimizer(GLPK.Optimizer(), Float64)

    lmin = zeros(T)
    lmax = 10.0 .* ones(T)
    price = collect(1:T)  # Cycle should start as early as possible

    # Instantiate initial model
    h = DR.DER.House(
        0, T,
        lmin, lmax, price,
        [l]
    )

    model, var2idx, con2idx = DR.DER.build_model!(h, model)

    MOI.optimize!(model)

    @test MOI.get(model, MOI.ObjectiveValue()) ≈ 8.0

    pnet = [
        var2idx[(:house, h.index, :pnet, t)]
        for t in 1:T
    ]

    @test MOI.get(model, MOI.VariablePrimal(), pnet[1]) ≈ 0.0
    @test MOI.get(model, MOI.VariablePrimal(), pnet[2]) ≈ 1.0
    @test MOI.get(model, MOI.VariablePrimal(), pnet[3]) ≈ 2.0
    @test MOI.get(model, MOI.VariablePrimal(), pnet[4]) ≈ 0.0


    return nothing
end

@testset "Shiftable" begin
    test_shiftableLoad()
end