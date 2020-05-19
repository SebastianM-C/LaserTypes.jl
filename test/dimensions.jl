using LaserTypes, Test
using Unitful
using StaticArrays
using UnitfulAtomic

@testset "$laser" for laser in [GaussLaser, LaguerreGaussLaser]
    units = [:SI_unitful, :atomic_unitful]
    x₀ = SVector{3}(0u"μm",0u"μm",0u"μm")
    t₀ = 0u"s"
    x = [x₀, auconvert.(x₀)]
    t = [t₀, auconvert(t₀)]

    @testset "$unit" for (unit, xᵢ, tᵢ) in zip(units, x, t)
        s = setup_laser(laser, unit)
        @test all(dimension.(E(xᵢ,tᵢ,s)) .== Ref(dimension(u"V/m")))
        @test all(dimension.(B(xᵢ,tᵢ,s)) .== Ref(dimension(u"T")))
    end
end
