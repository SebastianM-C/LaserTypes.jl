using LaserTypes, Test
using StaticArrays
using Unitful

@testset "Poynting vector" begin
    laser = setup_laser(GaussLaser, :SI_unitful, profile=ConstantProfile())
    r = SVector{3}(0,0,0)u"m"
    t = 0u"s"
    s = S(r, t, laser)

    @testset "Units" begin
        @test dimension(eltype(s)) == dimension(u"W/m^2")
    end
    @testset "Generic" begin
        @test isreal(s)
        @test ustrip(s[1]) isa Real
    end
end
