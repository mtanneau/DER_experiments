import GLPKMathProgInterface: GLPKSolverMIP

"""
    test_fixedLoad()

Run all unit tests for FixedLoad.
"""
function test_fixedLoad()

    @testset "FixedLoad" begin

        # Constructors
        test_fixedLoad_constructor()
        
        # Oracle
        test_fixedLoad_oracle()
    end

    return nothing
end

"""
    test_fixedLoad_constructor()

Run unit tests on constructors for FixedLoad.
"""
function test_fixedLoad_constructor()
    @testset "Constructors" begin

        # Constructor
        @test try
            DR.DER.FixedLoad(index=1, T=2, dt=1.0, load=[1.0, 2.0])
            true
        catch err
            false
        end

        @test try
            # Should raise a DimensionMismatch exception
            DR.DER.FixedLoad(T=1, load=[0.0, 1.0])
            false
        catch err
            isa(err, DimensionMismatch)
        end

        @test try
            # Should raise a DomainError
            # `dt` must be positive
            DR.DER.FixedLoad(index=0, T=1, dt=0.0, load=[1.0])
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
function test_fixedLoad_oracle()

    T = 2
    load = [1.0, 3.0]

    @testset "Oracle" begin
        l = DR.DER.FixedLoad(index=0, T=T, load=load)
        mip_solver = GLPKSolverMIP()

        # Instanciate oracle
        o = Linda.Oracle.LindaOracleMIP(l, mip_solver)
        
        # Query oracle with zero shadow prices
        Linda.Oracle.query!(o, zeros(2*T), 0.0)
        @test Linda.Oracle.get_sp_dual_bound(o) ≈ 0.0
        @test Linda.Oracle.get_num_new_columns(o) >= 1
        col = Linda.Oracle.get_new_columns(o)[1]
        @test col.col == [load; load]

        # Query oracle with non-zero shadow prices
        Linda.Oracle.query!(o, [ones(T); zeros(T)], 0.0)
        @test Linda.Oracle.get_sp_dual_bound(o) ≈ -sum(load)
        @test Linda.Oracle.get_num_new_columns(o) >= 1
        col = Linda.Oracle.get_new_columns(o)[1]
        @test col.col == [load; load]

    end

    return nothing
    
end

test_fixedLoad()