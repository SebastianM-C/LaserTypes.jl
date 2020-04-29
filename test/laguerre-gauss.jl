using LaserTypes, Test
using Unitful
using StaticArrays
using UnitfulAtomic
using Parameters
import PhysicalConstants.CODATA2018: c_0, m_e, e, α

Ex = LaguerreGauss.Ex
Ey = LaguerreGauss.Ey
Ez = LaguerreGauss.Ez
# Bx = LaguerreGauss.Bx
# By = LaguerreGauss.By
# Bz = LaguerreGauss.Bz

@testset "SI units" begin
    p = LaguerreGauss.LaguerreGaussParams()
    @unpack c, z_F, z_R, k, w₀, E₀ = p
    wz = LaserTypes.w(z_F, p)

    @test unit(c) == u"m/s"
    @test unit(k) == u"μm^-1"
    @test unit(z_R) == u"μm"

    x₀ = SVector{3}(0u"μm",0u"μm",z_F)
    t₀ = 0u"s"


    @testset "Dimensions" begin
        @test all(dimension.(E(x₀,t₀,p)) .== Ref(dimension(u"V/m")))
        # @test all(dimension.(B(x₀,t₀,p)) .== Ref(dimension(u"T")))
    end

    @testset "Units" begin
        @test all(unit.(E(x₀,t₀,p)) .== unit(p.E₀))
        # @test all(unit.(B(x₀,t₀,p)) .== unit(p.E₀/c))
    end

    @testset "Values at origin" begin
        @test_broken iszero(Ex(0u"μm", 0u"μm", 0u"μm", 0u"μm", p))
        @test iszero(Ey(0u"μm", 0u"μm", 0u"μm", 0u"μm", p))
        @test_broken iszero(Ez(0u"μm", 0u"μm", 0u"μm", 0u"μm", p))

        # @test iszero(Bx(0u"μm", 0u"μm", p))
        # @test By(0u"μm", 0u"μm", p) ≈ E₀ / c
        # @test iszero(Bz(0u"μm", 0u"μm", 0u"μm", 0u"μm", p))
    end
end
