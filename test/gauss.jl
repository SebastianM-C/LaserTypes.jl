using LaserTypes, Test
using LaserTypes: immutable_cache
using Unitful
using StaticArrays
using UnitfulAtomic

@testset "Values at origin" begin
    units = [:SI_unitful, :atomic_unitful, :SI, :atomic]
    x₀ = SVector{3}(0u"μm",0u"μm",0u"μm")
    t₀ = 0u"s"
    x = (x₀, auconvert.(x₀), ustrip.(x₀), ustrip.(auconvert.(x₀)))

    @testset "$unit" for (unit, xᵢ) in zip(units, x)
        s = setup_laser(GaussLaser, unit, profile=ConstantProfile())
        Ex, Ey, Ez = E(xᵢ, s)
        @test Ex ≈ immutable_cache(s, :E₀)
        @test iszero(Ey)
        @test iszero(Ez)

        Bx, By, Bz = B(xᵢ, s)
        @test iszero(Bx)
        @test By ≈ immutable_cache(s, :E₀) * immutable_cache(s, :inv_c)
        @test iszero(Bz)
    end
end

@testset "Orientation" begin
    s = setup_laser(GaussLaser, :SI, orientation=(:y,:x))
    x₀ = SVector{3}(0,0,0)
    t₀ = 0
    Ex, Ey, Ez = E(x₀, t₀, s)
    @test iszero(Ex)
    @test Ey ≈ immutable_cache(s, :E₀)
    @test iszero(Ez)

    Bx, By, Bz = B(x₀, t₀, s)
    @test iszero(Bx)
    @test iszero(By)
    @test Bz ≈ immutable_cache(s, :E₀) * immutable_cache(s, :inv_c)
end
