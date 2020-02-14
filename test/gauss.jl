using LaserTypes, Test
using Unitful
using StaticArrays
using UnitfulAtomic
import PhysicalConstants.CODATA2018: c_0, m_e, e, α

@testset "SI units" begin
    p = Gauss.LaserParams()
    @test dimension(p.c) == dimension(u"m/s")

    x₀ = SVector{3}(1,1,0.)u"m"

    E = Gauss.E
    B = Gauss.B
    @test all(dimension.(E(x₀,0u"s",p)) .== Ref(dimension(u"V/m")))
    @test all(dimension.(B(x₀,0u"s",p)) .== Ref(dimension(u"T")))
end

@testset "Atomic units" begin
    c = auconvert(c_0)
    q = auconvert(-e)
    m = auconvert(m_e)
    λ = auconvert(800u"nm")
    w0 = auconvert(58u"μm")
    τ0 = auconvert(18u"fs")
    p = Gauss.LaserParams(c=c, q=q, m_q=m, λ₀=λ, w₀=w0, τ₀=τ0)

    @test p.ω₀ ≈ 0.05695419u"Eh_au/ħ_au"

    x₀ = SVector{3}(1,1,0.)u"a0_au"

    E = Gauss.E
    B = Gauss.B
    @test all(dimension.(E(x₀,0u"ħ_au/Eh_au",p)) .== Ref(dimension(u"V/m")))
    @test all(dimension.(B(x₀,0u"ħ_au/Eh_au",p)) .== Ref(dimension(u"T")))
end

@testset "No units" begin
    c = 1 / α
    q = -1
    m = 1
    λ = austrip(800u"nm")
    w0 = austrip(58u"μm")
    τ0 = austrip(18u"fs")
    p = Gauss.LaserParams(c=c, q=q, m_q=m, λ₀=λ, w₀=w0, τ₀=τ0)

    @test p.ω₀ ≈ 0.05695419

    x₀ = SVector{3}(1,1,0.)

    E = Gauss.E
    B = Gauss.B
    @test E(x₀,0,p) ≈ [8.7831165e-7, 0, 0]
    @test B(x₀,0,p) ≈ [0, 6.409349797e-9, 0]
end
