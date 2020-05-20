using LaserTypes, Test
using Unitful
using StaticArrays
using UnitfulAtomic
using LinearAlgebra
import PhysicalConstants.CODATA2018: c_0

@testset "$laser" for laser in [GaussLaser, LaguerreGaussLaser]
    units = [:SI_unitful, :atomic_unitful, :SI, :atomic]
    x_si = SVector{4}(uconvert(u"μm", c_0*0u"s"), 0u"μm",0u"μm",1u"μm")
    x = (x_si, auconvert.(x_si), ustrip.(x_si), ustrip.(auconvert.(x_si)))

    @testset "$unit" for (unit, xᵢ) in zip(units, x)
        s = setup_laser(laser, unit, profile=ConstantProfile())
        F = Fμν(xᵢ,s)
        c = s.c

        if unit ∈ [:SI_unitful, :atomic_unitful]
            @test all(dimension.(F) .== Ref(dimension(u"T")))
        end
        @test size(F) == (4, 4)
        @test det(F) ≈ 1/c^2 * (E(xᵢ[2:4], xᵢ[1]/c, s) ⋅ B(xᵢ[2:4], xᵢ[1]/c, s))^2
    end
end
