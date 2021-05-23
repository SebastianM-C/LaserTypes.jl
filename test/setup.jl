using LaserTypes, Test
using Unitful
using UnitfulAtomic
import PhysicalConstants.CODATA2018: c_0, m_e, μ_0, e, α

@testset "SI Unitful" begin
    c = c_0
    q = -e
    m = m_e
    λ = 0.8u"μm"
    w0 = 58.0u"μm"
    τ = 18.0u"fs"
    a₀ = 2.0
    p = GaussLaser(:SI_unitful; λ, w₀=w0, a₀=2.0, profile=GaussProfile(;τ,z₀=0.0u"μm"))
    s = setup_laser(GaussLaser, :SI_unitful; a₀=2.0, profile = GaussProfile, τ)
    @test p == s

    @test unit(s.constants.c) == u"m/s"
    @test unit(s.derived.k) == u"μm^-1"
    @test unit(s.derived.z_R) == u"μm"
end

@testset "Atomic Unitful" begin
    c = auconvert(c_0)
    q = auconvert(-e)
    m = auconvert(m_e)
    μ₀ = auconvert(μ_0)
    λ = auconvert(800u"nm")
    w0 = auconvert(58u"μm")
    τ = auconvert(18u"fs")
    a₀=2.0

    p = GaussLaser(:atomic_unitful; λ, w₀=w0, a₀=2.0, profile=GaussProfile(;τ,z₀=0.0u"a0_au"))
    s = setup_laser(GaussLaser, :atomic_unitful; a₀=2.0, profile = GaussProfile, τ)
    @test p == s

    @test unit(s.constants.c) == u"a0_au*Eh_au/ħ_au"
    @test unit(s.derived.k) == u"a0_au^-1"
    @test unit(s.derived.z_R) == u"a0_au"
    @test s.derived.ω ≈ 0.05695419u"Eh_au/ħ_au"
end

@testset "SI" begin
    c = ustrip(c_0)
    q = ustrip(-e)
    m = ustrip(m_e)
    μ₀ = ustrip(μ_0)
    λ = ustrip(u"m", 0.8u"μm")
    w0 = ustrip(u"m", 58.0u"μm")
    τ = ustrip(u"s", 18.0u"fs")
    p = GaussLaser(:SI; λ, w₀=w0, a₀=2.0, profile=GaussProfile(;τ, z₀=0.0))
    s = setup_laser(GaussLaser, :SI, a₀=2.0)
    @test p == s

    @test s.derived.ω ≈ 2.354564459e15
end

@testset "Atomic" begin
    c = ustrip(auconvert(c_0))
    q = ustrip(auconvert(-e))
    m = ustrip(auconvert(m_e))
    μ₀ = ustrip(auconvert(μ_0))
    λ = ustrip(auconvert(800u"nm"))
    w0 = ustrip(auconvert(58u"μm"))
    τ = ustrip(auconvert(18u"fs"))
    a₀=2.0

    p = GaussLaser(:atomic; λ, w₀=w0, a₀=2.0, profile=GaussProfile(τ=τ, z₀=0.0))
    s = setup_laser(GaussLaser, :atomic, a₀=2.0)
    @test p == s

    @test s.derived.ω ≈ 0.05695419
end
