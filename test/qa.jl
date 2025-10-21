using LaserTypes, Test
using JET

@testset "JET type stability tests" begin
    @testset "$laser" for laser in [GaussLaser, LaguerreGaussLaser]
        # Setup laser with atomic units (simplest for testing)
        if laser == GaussLaser
            lt_laser = setup_laser(GaussLaser, :atomic, profile=ConstantProfile())
        else
            lt_laser = setup_laser(LaguerreGaussLaser, :atomic, profile=ConstantProfile(), p=1, m=1)
        end

        r = [0.0, 0.0, 0.0]
        t = 0.0

        @testset "E field - Val(:real)" begin
            @test_opt LaserTypes.E(r, t, lt_laser, Val(:real))
        end

        @testset "E field - Val(:complex)" begin
            @test_opt LaserTypes.E(r, t, lt_laser, Val(:complex))
        end

        @testset "B field - Val(:real)" begin
            @test_opt LaserTypes.B(r, t, lt_laser, Val(:real))
        end

        @testset "B field - Val(:complex)" begin
            @test_opt LaserTypes.B(r, t, lt_laser, Val(:complex))
        end

        @testset "EB field - Val(:real)" begin
            @test_opt LaserTypes.EB(r, t, lt_laser, Val(:real))
        end

        @testset "EB field - Val(:complex)" begin
            @test_opt LaserTypes.EB(r, t, lt_laser, Val(:complex))
        end
    end
end
