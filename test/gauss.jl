using LaserTypes, Test
using Unitful
using StaticArrays
using UnitfulAtomic

@testset "Values at origin" begin
    units = [:SI_unitful, :atomic_unitful, :SI, :atomic]
    x₀ = SVector{3}(0u"μm",0u"μm",0u"μm")
    t₀ = 0u"s"
    x = (x₀, auconvert.(x₀), ustrip.(x₀), ustrip.(auconvert.(x₀)))
    t = (t₀, auconvert(t₀), ustrip(t₀), ustrip(auconvert(t₀)))

    @testset "$unit" for (unit, xᵢ, tᵢ) in zip(units, x, t)
        s = setup_laser(GaussLaser, unit, profile=ConstantProfile())
        Ex, Ey, Ez = E(xᵢ, tᵢ, s)
        @test Ex ≈ s.E₀
        @test iszero(Ey)
        @test iszero(Ez)

        Bx, By, Bz = B(xᵢ, tᵢ, s)
        @test iszero(Bx)
        @test By ≈ s.E₀ / s.c
        @test iszero(Bz)
    end
end
