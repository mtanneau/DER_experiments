"""
    test_curtailableLoad()

Run all unit tests for CurtailableLoad.
"""
function test_curtailableLoad()

    @testset "Curt." begin
        # Constructor
        test_curtailableLoad_constructor()

        # Oracle
        test_curtailableLoad_oracle()

    end
    return nothing
end

"""
    test_curtailableLoad_constructor()

Run unit tests on constructors for FixedLoad.
"""
function test_curtailableLoad_constructor()
    @testset "Constructors" begin

        # Constructor
        @test try
            DR.DER.CurtailableLoad(
                index=1, T=2, dt=1.0,
                load=[1.0, 2.0], binaryFlag=false
            )
            true
        catch err
            false
        end

        @test try
            # Should raise a DimensionMismatch exception
            DR.DER.CurtailableLoad(T=1, load=[0.0, 1.0])
            false
        catch err
            isa(err, DimensionMismatch)
        end

        @test try
            # Should raise a DomainError
            # `dt` must be positive
            DR.DER.CurtailableLoad(index=0, T=1, dt=0.0, load=[1.0])
            false
        catch err
            isa(err, DomainError)
        end

    end
    return nothing
end

"""
    test_fixedLoad_oracle()

Run unit tests on oracle instanciation for FixedLoad
"""
function test_curtailableLoad_oracle()

    T = 2
    load = [1.0, 3.0]

    @testset "Oracle" begin
        l = DR.DER.CurtailableLoad(index=0, T=T, load=load)
        mip_solver = GLPKSolverMIP()

        # Instanciate oracle
        o = Linda.Oracle.LindaOracleMIP(l, mip_solver)
        
        # Query oracle with zero shadow prices
        Linda.Oracle.query!(o, zeros(2*T), 0.0)
        @test Linda.Oracle.get_sp_dual_bound(o) ≈ 0.0
        @test Linda.Oracle.get_num_new_columns(o) >= 1

        # Query oracle with positive shadow prices
        # (we're minimizing cost so consumption is maximized)
        Linda.Oracle.query!(o, [ones(T); zeros(T)], 0.0)
        @test Linda.Oracle.get_sp_dual_bound(o) ≈ -sum(load)
        @test Linda.Oracle.get_num_new_columns(o) >= 1
        col = Linda.Oracle.get_new_columns(o)[1]
        @test col.col == [load; load]

        # Query oracle with negative shadow prices
        # (we're minimizing cost so consumption is minimized)
        Linda.Oracle.query!(o, [-ones(T); zeros(T)], 0.0)
        @test Linda.Oracle.get_sp_dual_bound(o) ≈ 0.0
        @test Linda.Oracle.get_num_new_columns(o) >= 1
        col = Linda.Oracle.get_new_columns(o)[1]
        @test col.col == [zeros(T); zeros(T)]

    end

    return nothing
    
end

test_curtailableLoad()