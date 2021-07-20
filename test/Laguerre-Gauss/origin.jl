using LaserTypes
using LaserTypes: pochhammer, immutable_cache
using Unitful, UnitfulAtomic
using StaticArrays
using UnPack
using Test

units = [:SI_unitful, :atomic_unitful, :SI, :atomic]
x₀ = SVector{3}(0u"μm",0u"μm",0u"μm")
t₀ = 0u"s"
x = (x₀, auconvert.(x₀), ustrip.(x₀), ustrip.(auconvert.(x₀)))
t = (t₀, auconvert.(t₀), ustrip.(t₀), ustrip.(auconvert.(t₀)))

@testset "$unit" for (unit, xᵢ, tᵢ) in zip(units, x, t)
    s = setup_laser(LaguerreGaussLaser, unit, m = 1, p = 1, profile=ConstantProfile)
    Ex, Ey, Ez = E(xᵢ, tᵢ, s)
    @test iszero(Ex)
    @test iszero(Ey)
    @test iszero(Ez)

    @unpack E₀, Nₚₘ, z_R = immutable_cache(s)
    c = s.constants.c
    Bx, By, Bz = B(xᵢ, tᵢ,s)
    @test iszero(Bx)
    @test iszero(By)
    @test Bz ≈ - factorial(s.p)/pochhammer(abs(s.m)+1,s.p)*E₀*Nₚₘ*(√2*s.w₀)/(c*z_R)
end

@testset "Default laser" begin
    s = setup_laser(LaguerreGaussLaser, :atomic; λ = 1e5, a₀ = 1, w₀ = 1e6,
        profile = GaussProfile(τ = 100, z₀ = 0.))

    @test E([0, 0, 0], 0, s)[1] == s.derived.E₀
end
