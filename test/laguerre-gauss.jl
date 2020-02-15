using LaserTypes, Test
using Unitful
using StaticArrays
using UnitfulAtomic
using Parameters
import PhysicalConstants.CODATA2018: c_0, m_e, e, α

@testset "SI units" begin
    p = Gauss.LaserParams()
    @unpack c, z_F, z_R, k, w₀, E₀ = p
    wz = LaserTypes.w(z_F, p)

    @test unit(c) == u"m/s"
    @test unit(k) == u"μm^-1"
    @test unit(z_R) == u"μm"
end
