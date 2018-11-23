"""
    test_thermalLoad()

Run all unit tests for CurtailableLoad.
"""
function test_thermalLoad()

    @testset "Thermal" begin
        # Constructor
        test_thermalLoad_constructor()

        # Oracle
        test_thermalLoad_oracle()

    end
    return nothing
end

"""
    test_thermalLoad_constructor()

Run unit tests on constructors for ThermalLoad.
"""
function test_thermalLoad_constructor()
    @testset "Constructors" begin

        # Constructor
        @test try
            DR.DER.ThermalLoad(
                index=1, T=2, dt=1.0,
                temp_min=[18.0, 18.0],
                temp_max=[22.0, 22.0],
                temp_ext=[0.0, 0.0],
                temp_init=20.0,
                pwr_min=[2.0, 2.0],
                pwr_max=[10.0, 10.0],
                C=1.0, η=1.0, μ=0.2,
                binary_flag=true
            )
            true
        catch err
            false
        end

        @test try
            # Should raise a DimensionMismatch exception
            DR.DER.ThermalLoad(temp_min=[18.0])
            DR.DER.ThermalLoad(temp_max=[18.0])
            DR.DER.ThermalLoad(temp_ext=[18.0])
            DR.DER.ThermalLoad(pwr_min=[18.0])
            DR.DER.ThermalLoad(pwr_max=[18.0])
            false
        catch err
            isa(err, DimensionMismatch)
        end

        @test try
            # Should raise a DomainError
            # `dt` must be positive
            DR.DER.ThermalLoad(dt=0.0)
            DR.DER.ThermalLoad(C=0.0)
            false
        catch err
            isa(err, DomainError)
        end

    end
    return nothing
end

"""
    test_thermalLoad_oracle()

Run unit tests on oracle instanciation for FixedLoad
"""
function test_thermalLoad_oracle()

    T = 2

    @testset "Oracle" begin
        l = DR.DER.ThermalLoad(
            index=1, T=2, dt=1.0,
            temp_min=[20.0, 20.0],
            temp_max=[20.0, 20.0],
            temp_ext=[0.0, 0.0],
            temp_init=20.0,
            pwr_min=[0.0, 0.0],
            pwr_max=[10.0, 10.0],
            C=1.0, η=1.0, μ=0.2,
            binary_flag=false
        )
        mip_solver = GLPKSolverMIP()

        # Instanciate oracle
        o = Linda.Oracle.LindaOracleMIP(l, mip_solver)
        
        # Query oracle with zero shadow prices
        Linda.Oracle.query!(o, zeros(2*T), 0.0)
        @test Linda.Oracle.get_sp_dual_bound(o) ≈ 0.0
        @test Linda.Oracle.get_num_new_columns(o) >= 1
        col = Linda.Oracle.get_new_columns(o)[1]
        @test col.col == [4.0, 4.0, 4.0, 4.0]

        # Query oracle with positive shadow prices
        # (we're minimizing cost so consumption is maximized)
        Linda.Oracle.query!(o, [ones(T); zeros(T)], 0.0)
        # @test Linda.Oracle.get_sp_dual_bound(o) ≈ -sum(load)
        @test Linda.Oracle.get_num_new_columns(o) >= 1
        

        # Query oracle with negative shadow prices
        # (we're minimizing cost so consumption is minimized)
        Linda.Oracle.query!(o, [-ones(T); zeros(T)], 0.0)
        # @test Linda.Oracle.get_sp_dual_bound(o) ≈ 0.0
        @test Linda.Oracle.get_num_new_columns(o) >= 1
        # col = Linda.Oracle.get_new_columns(o)[1]
        # @test col.col == [zeros(T); zeros(T)]

    end

    return nothing
    
end

test_thermalLoad()