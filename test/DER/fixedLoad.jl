function test_fixedLoad()

    T = 2
    load = [1.0, 3.0]

    l = DR.DER.FixedLoad(index=0, T=T, load=load)
    model = GLPK.Optimizer()

    # Instantiate initial model
    h = DR.DER.House(
        0, T,
        [0.0, 0.0], [10.0, 10.0],
        [1.0, 1.0],
        [l]
    )

    model, var2idx, con2idx = DR.DER.build_model!(h, model)

    MOI.optimize!(model)

    # Test solution
    @test MOI.get(model, MOI.ObjectiveValue()) ≈ 4.0
    for t in 1:T
        pnet = var2idx[(:house, h.index, :pnet, t)]
        pfix = var2idx[(:fixed, l.index, :pnet, t)]
        @test MOI.get(model, MOI.VariablePrimal(), pnet) ≈ load[t]
        @test MOI.get(model, MOI.VariablePrimal(), pfix) ≈ load[t]
    end

    return nothing
end

@testset "FixedLoad" begin
    test_fixedLoad()
end