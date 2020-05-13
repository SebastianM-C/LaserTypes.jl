using LaserTypes, Test
using Unitful
using StaticArrays
using UnitfulAtomic
using Parameters
import PhysicalConstants.CODATA2018: c_0, m_e, e, α

@testset "SI units" begin
    p = LaguerreGaussLaser()
    @unpack c, profile, z_R, k, w₀, E₀ = p
    z₀ = profile.z₀
    wz = LaserTypes.w(z₀, p)

    @test unit(c) == u"m/s"
    @test unit(k) == u"μm^-1"
    @test unit(z_R) == u"μm"

    x₀ = SVector{3}(0u"μm",0u"μm",z₀)
    t₀ = 0u"s"

    @testset "Dimensions" begin
        @test all(dimension.(E(x₀,t₀,p)) .== Ref(dimension(u"V/m")))
        @test all(dimension.(B(x₀,t₀,p)) .== Ref(dimension(u"T")))
    end

    @testset "Values at origin" begin
        Ex, Ey, Ez = E(0u"μm", 0u"μm", 0u"μm", p)
        @test Ex ≈ E₀
        @test iszero(Ey)
        @test iszero(Ez)

        Bx, By, Bz = B(0u"μm", 0u"μm", 0u"μm", p)
        @test iszero(Bx)
        @test By ≈ E₀ / c
        @test iszero(Bz)
    end
end
