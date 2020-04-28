using LaserTypes, Test
using Unitful
using StaticArrays
using UnitfulAtomic
using Parameters
import PhysicalConstants.CODATA2018: c_0, m_e, e, α

E = Gauss.E
B = Gauss.B

Ex = Gauss.Ex
Ey = Gauss.Ey
Ez = Gauss.Ez
Bx = Gauss.Bx
By = Gauss.By
Bz = Gauss.Bz

@testset "SI units" begin
    p = Gauss.LaserParams()
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
        @test Ex(0u"μm", 0u"μm", p) ≈ E₀
        @test iszero(Ey(0u"μm", 0u"μm", p))
        @test iszero(Ez(0u"μm", 0u"μm", 0u"μm", 0u"μm", p))

        @test iszero(Bx(0u"μm", 0u"μm", p))
        @test By(0u"μm", 0u"μm", p) ≈ E₀ / c
        @test iszero(Bz(0u"μm", 0u"μm", 0u"μm", 0u"μm", p))
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

    p = Gauss.LaserParams(c=c, q=q, m_q=m, λ=λ, w₀=w0, τ₀=τ0)
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
        @test Ex(0u"a0_au", 0u"a0_au", p) ≈ E₀
        @test iszero(Ey(0u"a0_au", 0u"a0_au", p))
        @test iszero(Ez(0u"a0_au", 0u"a0_au", 0u"a0_au", 0u"a0_au", p))

        @test iszero(Bx(0u"a0_au", 0u"a0_au", p))
        @test By(0u"a0_au", 0u"a0_au", p) ≈ E₀ / c
        @test iszero(Bz(0u"a0_au", 0u"a0_au", 0u"a0_au", 0u"a0_au", p))
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

    p = Gauss.LaserParams(c=c, q=q, m_q=m, λ=λ, w₀=w0, τ₀=τ0)
    @unpack c, z_F, z_R, k, w₀, E₀ = p
    wz = LaserTypes.w(z_F, p)

    @test p.ω ≈ 0.05695419

    x₀ = SVector{3}(0, 0, z_F)
    t₀ = 0.

    E = Gauss.E
    B = Gauss.B

    @testset "Values at origin" begin
        @test Ex(0, 0, p) ≈ E₀
        @test iszero(Ey(0, 0, p))
        @test iszero(Ez(0, 0, 0, 0, p))

        @test iszero(Bx(0, 0, p))
        @test By(0, 0, p) ≈ E₀ / c
        @test iszero(Bz(0, 0, 0, 0, p))
    end

    @testset "Values at z_F" begin
        @test E(x₀,t₀,p) ≈ [E₀*w₀/wz*real(exp(-im*k*z_F+im*atan(z_F,z_R))), 0, 0]
        @test B(x₀,t₀,p) ≈ [0, E₀*w₀/wz*real(exp(-im*k*z_F+im*atan(z_F,z_R)))/c, 0]
    end
end
