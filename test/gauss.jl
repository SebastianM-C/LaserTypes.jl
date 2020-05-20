using LaserTypes, Test
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
        @test Ex ≈ s.E₀
        @test iszero(Ey)
        @test iszero(Ez)

        Bx, By, Bz = B(xᵢ, s)
        @test iszero(Bx)
        @test By ≈ s.E₀ / s.c
        @test iszero(Bz)
    end
end
