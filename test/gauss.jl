using LaserTypes, Test
using Unitful
using StaticArrays
using UnitfulAtomic
using Parameters
import PhysicalConstants.CODATA2018: c_0, m_e, e, α

@testset "SI units" begin
    p = GaussLaser()
    @unpack c, z_F, z_R, k, w₀, E₀ = p
    wz = LaserTypes.w(z_F, p)

    @test unit(c) == u"m/s"
    @test unit(k) == u"μm^-1"
    @test unit(z_R) == u"μm"

    x₀ = SVector{3}(0u"μm",0u"μm",z_F)
    t₀ = 0u"s"

    @testset "Dimensions" begin
        @test all(dimension.(E(x₀,t₀,p)) .== Ref(dimension(u"V/m")))
        @test all(dimension.(B(x₀,t₀,p)) .== Ref(dimension(u"T")))
    end

    @testset "Units" begin
        @test all(unit.(E(x₀,t₀,p)) .== unit(p.E₀))
        @test all(unit.(B(x₀,t₀,p)) .== unit(p.E₀/c))
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

    @testset "Values at z_F" begin
        @test E(x₀,t₀,p) ≈ [E₀*w₀/wz*real(exp(-im*k*z_F+im*atan(z_F,z_R))), 0u"V/m", 0u"V/m"]
        @test B(x₀,t₀,p) ≈ [0u"T", E₀*w₀/wz*real(exp(-im*k*z_F+im*atan(z_F,z_R)))/c, 0u"T"]
    end
end

@testset "Atomic units" begin
    c = auconvert(c_0)
    q = auconvert(-e)
    m = auconvert(m_e)
    λ = auconvert(800u"nm")
    w0 = auconvert(58u"μm")
    τ0 = auconvert(18u"fs")

    p = GaussLaser(c=c, q=q, m_q=m, λ=λ, w₀=w0, τ₀=τ0)
    @unpack c, z_F, z_R, k, w₀, E₀ = p
    wz = LaserTypes.w(z_F, p)

    @test unit(c) == u"a0_au*Eh_au/ħ_au"
    @test unit(k) == u"a0_au^-1"
    @test unit(z_R) == u"a0_au"

    @test p.ω ≈ 0.05695419u"Eh_au/ħ_au"

    x₀ = SVector{3}(0u"a0_au",0u"a0_au",z_F)
    t₀ = 0u"ħ_au/Eh_au"

    @testset "Dimensions" begin
        @test all(dimension.(E(x₀,t₀,p)) .== Ref(dimension(u"V/m")))
        @test all(dimension.(B(x₀,t₀,p)) .== Ref(dimension(u"T")))
    end

    @testset "Units" begin
        @test all(unit.(E(x₀,t₀,p)) .== unit(p.E₀))
        @test all(unit.(B(x₀,t₀,p)) .== unit(p.E₀/c))
    end

    @testset "Values at origin" begin
        Ex, Ey, Ez = E(0u"a0_au", 0u"a0_au", 0u"a0_au", p)
        @test Ex ≈ E₀
        @test iszero(Ey)
        @test iszero(Ez)

        Bx, By, Bz = B(0u"a0_au", 0u"a0_au", 0u"a0_au", p)
        @test iszero(Bx)
        @test By ≈ E₀ / c
        @test iszero(Bz)
    end

    @testset "Values at z_F" begin
        @test E(x₀,t₀,p) ≈ [E₀*w₀/wz*real(exp(-im*k*z_F+im*atan(z_F,z_R))), 0u"V/m", 0u"V/m"]
        @test B(x₀,t₀,p) ≈ [0u"T", E₀*w₀/wz*real(exp(-im*k*z_F+im*atan(z_F,z_R)))/c, 0u"T"]
    end
end

@testset "No units" begin
    c = 1 / α
    q = -1
    m = 1
    λ = austrip(800u"nm")
    w0 = austrip(58u"μm")
    τ0 = austrip(18u"fs")

    p = GaussLaser(c=c, q=q, m_q=m, λ=λ, w₀=w0, τ₀=τ0)
    @unpack c, z_F, z_R, k, w₀, E₀ = p
    wz = LaserTypes.w(z_F, p)

    @test p.ω ≈ 0.05695419

    x₀ = SVector{3}(0, 0, z_F)
    t₀ = 0.

    @testset "Values at origin" begin
        Ex, Ey, Ez = E(0, 0, 0, p)
        @test Ex ≈ E₀
        @test iszero(Ey)
        @test iszero(Ez)

        Bx, By, Bz = B(0, 0, 0, p)
        @test iszero(Bx)
        @test By ≈ E₀ / c
        @test iszero(Bz)
    end

    @testset "Values at z_F" begin
        @test E(x₀,t₀,p) ≈ [E₀*w₀/wz*real(exp(-im*k*z_F+im*atan(z_F,z_R))), 0, 0]
        @test B(x₀,t₀,p) ≈ [0, E₀*w₀/wz*real(exp(-im*k*z_F+im*atan(z_F,z_R)))/c, 0]
    end
end
