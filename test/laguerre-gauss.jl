using LaserTypes, Test
using LaserTypes: pochhammer, immutable_cache
using UnPack
using Unitful
using StaticArrays
using UnitfulAtomic
using LinearAlgebra


@testset "Values at origin" begin
    units = [:SI_unitful, :atomic_unitful, :SI, :atomic]
    x₀ = SVector{3}(0u"μm",0u"μm",0u"μm")
    t₀ = 0u"s"
    x = (x₀, auconvert.(x₀), ustrip.(x₀), ustrip.(auconvert.(x₀)))
    t = (t₀, auconvert.(t₀), ustrip.(t₀), ustrip.(auconvert.(t₀)))

    @testset "$unit" for (unit, xᵢ, tᵢ) in zip(units, x, t)
        s = setup_laser(LaguerreGaussLaser, unit, m = 1, p = 1, profile=ConstantProfile())
        Ex, Ey, Ez = E(xᵢ, tᵢ, s)
        @test iszero(Ex)
        @test iszero(Ey)
        @test iszero(Ez)

        @unpack E₀, Nₚₘ, z_R = immutable_cache(s)
        c = s.constants.c
        Bx, By, Bz = B(xᵢ, tᵢ,s)
        @test iszero(Bx)
        @test iszero(By)
        @test isapprox(Bz, - factorial(s.p)/pochhammer(abs(s.m)+1,s.p)*E₀*Nₚₘ*(√2*s.w₀)/(c*z_R))
    end

    @testset "Default laser" begin
        s = setup_laser(LaguerreGaussLaser, :atomic; λ = 1e5, a₀ = 1, w₀ = 1e6,
            profile = GaussProfile(τ = 100, z₀ = 0.))

        @test E([0, 0, 0], 0, s)[1] == s.derived.E₀
    end
end

@testset "Divergence in Polynomial Roots" begin
    ω = 0.057
    T₀ = 2π/ω
    c = 137.036
    λ = c*T₀
    w₀ = 75 * λ
    a₀ = 2.
    s = setup_laser(LaguerreGaussLaser, :atomic, m = 1, p = 1, λ = λ, a₀ = a₀, w₀ = w₀,
        profile=ConstantProfile())

    fieldE(x,y,z,t) = E([x,y,z],t,s)
    fieldB(x,y,z,t) = B([x,y,z],t,s)

    wenergy(x,y,t) = fieldE(x,y,0,t)⋅fieldE(x,y,0,t) + c*c * fieldB(x,y,0,t)⋅fieldB(x,y,0,t)
    wenergy₀(x,y) = wenergy(x,y,0)

    domainXY = -4*w₀:λ:4*w₀
    wXY = map(r -> wenergy₀(r[1],r[2]), Iterators.product(domainXY,domainXY))
    wMax = maximum(wXY)
    @test wMax < 1e6
end

@testset "LG(p = 0, m = 0, x, y, z, t) == G(x, y, z, t)" begin
    ω = 0.057
    T₀ = 2π/ω
    c = 137.036
    λ = c*T₀
    w₀ = 75 * λ
    a₀ = 2.
    LG = setup_laser(LaguerreGaussLaser, :atomic, m = 0, p = 0, λ = λ, a₀ = a₀, w₀ = w₀,
        profile=ConstantProfile())
    G = setup_laser(GaussLaser, :atomic, λ = λ, a₀ = a₀, w₀ = w₀,
        profile=ConstantProfile())
    for i in 1:100
        x, y, z = w₀ * (0.5 .- rand(3))
        t = T₀ * rand()
        @test Fμν([t*c, x, y, z], G) ≈ Fμν([t*c, x, y, z], LG) rtol = 1e-8
    end
end
