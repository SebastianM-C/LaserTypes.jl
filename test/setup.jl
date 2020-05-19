using LaserTypes, Test
using Unitful
using UnitfulAtomic
import PhysicalConstants.CODATA2018: c_0, m_e, e, α

@testset "SI Unitful" begin
    c = c_0
    q = -e
    m = m_e
    λ = 0.8u"μm"
    w0 = 58.0u"μm"
    τ = 18.0u"fs"
    p = GaussLaser(c=c, q=q, m_q=m, λ=λ, w₀=w0, profile=GaussProfile(c=c,τ=τ))
    s = setup_laser(GaussLaser, :SI_unitful)
    @test p == s

    @test unit(s.c) == u"m/s"
    @test unit(s.k) == u"μm^-1"
    @test unit(s.z_R) == u"μm"
end

@testset "Atomic Unitful" begin
    c = auconvert(c_0)
    q = auconvert(-e)
    m = auconvert(m_e)
    λ = auconvert(800u"nm")
    w0 = auconvert(58u"μm")
    τ = auconvert(18u"fs")

    p = GaussLaser(c=c, q=q, m_q=m, λ=λ, w₀=w0, profile=GaussProfile(c=c,τ=τ))
    s = setup_laser(GaussLaser, :atomic_unitful)
    @test p == s

    @test unit(s.c) == u"a0_au*Eh_au/ħ_au"
    @test unit(s.k) == u"a0_au^-1"
    @test unit(s.z_R) == u"a0_au"
    @test s.ω ≈ 0.05695419u"Eh_au/ħ_au"
end

@testset "SI" begin
    c = ustrip(c_0)
    q = ustrip(-e)
    m = ustrip(m_e)
    λ = ustrip(u"m", 0.8u"μm")
    w0 = ustrip(u"m", 58.0u"μm")
    τ = ustrip(u"s", 18.0u"fs")
    p = GaussLaser(c=c, q=q, m_q=m, λ=λ, w₀=w0, profile=GaussProfile(c=c,τ=τ))
    s = setup_laser(GaussLaser, :SI)
    @test p == s

    @test s.ω ≈ 2.354564459e15
end

@testset "Atomic" begin
    c = ustrip(auconvert(c_0))
    q = ustrip(auconvert(-e))
    m = ustrip(auconvert(m_e))
    λ = ustrip(auconvert(800u"nm"))
    w0 = ustrip(auconvert(58u"μm"))
    τ = ustrip(auconvert(18u"fs"))

    p = GaussLaser(c=c, q=q, m_q=m, λ=λ, w₀=w0, profile=GaussProfile(c=c,τ=τ))
    s = setup_laser(GaussLaser, :atomic)
    @test p == s

    @test s.ω ≈ 0.05695419
end