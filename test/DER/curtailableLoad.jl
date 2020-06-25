function test_curtailableLoad()

    T = 2
    load = [2.0, 3.0]

    l = DR.DER.CurtailableLoad(index=0, T=T, load=load, binflag=false)

    model = GLPK.Optimizer()

    lmin = [1.0, 0.0]
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

    @test MOI.get(model, MOI.ObjectiveValue()) ≈ 1.0-0.03

    for t in 1:T
        pnet = var2idx[(:house, h.index, :pnet, t)]
        pres = var2idx[(:curt, l.index, :pnet, t)]
        if price[t] >= 0.0
            @test MOI.get(model, MOI.VariablePrimal(), pnet) ≈ lmin[t]
            @test MOI.get(model, MOI.VariablePrimal(), pres) ≈ lmin[t]
        else
            @test MOI.get(model, MOI.VariablePrimal(), pnet) ≈ min(lmax[t], load[t])
            @test MOI.get(model, MOI.VariablePrimal(), pres) ≈ min(lmax[t], load[t])
        end
    end

    return nothing
    
end

@testset "Curt." begin
    test_curtailableLoad()
end