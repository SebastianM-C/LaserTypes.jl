using LaserTypes
using LaserTypes: pochhammer, immutable_cache
using StaticArrays
using UnPack
using Test

x₀ = SVector{3}(0,0,0)
t₀ = 0

@testset "(:y, :x) orientation" begin
    s = setup_laser(LaguerreGaussLaser, :SI, m = 1, p = 1, profile=ConstantProfile, orientation=(:y,:x))
    Ex, Ey, Ez = E(x₀, t₀, s)
    @test iszero(Ex)
    @test iszero(Ey)
    @test iszero(Ez)

    @unpack E₀, Nₚₘ, z_R = immutable_cache(s)
    c = s.constants.c
    Bx, By, Bz = B(x₀, t₀, s)
    @test Bx ≈ - factorial(s.p)/pochhammer(abs(s.m)+1,s.p)*E₀*Nₚₘ*(√2*s.w₀)/(c*z_R)
    @test iszero(By)
    @test iszero(Bz)
end
