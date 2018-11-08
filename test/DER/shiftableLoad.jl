"""
    test_shiftableLoad()

Run all unit tests for ShiftableLoad.
"""
function test_shiftableLoad()

    @testset "Shift." begin
        # Constructor
        test_shiftableLoad_constructor()

        # Oracle
        test_shiftableLoad_oracle()

    end
    return nothing
end

"""
    test_shiftableLoad_constructor()

Run unit tests on constructors for ShiftableLoad.
"""
function test_shiftableLoad_constructor()
    @testset "Constructors" begin

        # Constructor
        @test try
            DR.DER.ShiftableLoad(
                T=2,
                cycle_begin_min=1,
                cycle_begin_max=2,
                cycle_length=1,
                cycle_load=[1.0]
            )
            true
        catch err
            rethrow(err)
            false
        end

        @test try
            # Should raise a DimensionMismatch exception
            DR.DER.ShiftableLoad(
                T=3,
                cycle_begin_min=1,
                cycle_begin_max=2,
                cycle_length=1,
                cycle_load=[1.0, 2.0]
            )
            false
        catch err
            isa(err, DimensionMismatch)
        end

    end
    return nothing
end

"""
    test_fixedLoad_oracle()

Run unit tests on oracle instanciation for ShiftableLoad
"""
function test_shiftableLoad_oracle()

    T = 4

    @testset "Oracle" begin
        l = DR.DER.ShiftableLoad(
            T=T,
            cycle_begin_min=1,
            cycle_begin_max=3,
            cycle_length=2,
            cycle_load=[1.0, 2.0]
        )
        mip_solver = GLPKSolverMIP()

        # Instanciate oracle
        o = Linda.Oracle.LindaOracleMIP(l, mip_solver)

        # Query oracle with zero shadow prices
        Linda.Oracle.query!(o, zeros(2*T), 0.0)
        @test Linda.Oracle.get_sp_dual_bound(o) ≈ 0.0
        @test Linda.Oracle.get_num_new_columns(o) >= 1


        # Query oracle with increasing shadow prices
        # (we're minimizing cost so consumption is maximized)
        Linda.Oracle.query!(o, [zeros(T); collect(1:T)], 0.0)
        # Cycle should be shifted at latest possible time
        # Total cost should be `- (1.0 * 3.0 + 2.0 * 4.0) = -11.0`
        @test Linda.Oracle.get_sp_dual_bound(o) ≈ -11.0
        @test Linda.Oracle.get_num_new_columns(o) >= 1
        col = Linda.Oracle.get_new_columns(o)[1]
        @test col.col == [[0.0, 0.0, 1.0, 2.0]; [0.0, 0.0, 1.0, 2.0]]

        # Query oracle with negative shadow prices
        # (we're minimizing cost so consumption is minimized)
        Linda.Oracle.query!(o, [-ones(T); zeros(T)], 0.0)
        @test Linda.Oracle.get_sp_dual_bound(o) ≈ sum(l.cycle_load)
        @test Linda.Oracle.get_num_new_columns(o) >= 1
        col = Linda.Oracle.get_new_columns(o)[1]

    end

    return nothing
    
end

test_shiftableLoad()