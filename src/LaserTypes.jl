# # LaserTypes

module LaserTypes

export E, B, GaussLaser, LaguerreGaussLaser, ConstantProfile, GaussProfile,
    QuasiRectangularProfile, setup_laser

using Parameters

using Unitful
using UnitfulAtomic
using Parameters
using GeometryTypes: Vec3
import PhysicalConstants.CODATA2018: c_0, m_e, e
using HypergeometricFunctions

const _₁F₁ = HypergeometricFunctions.drummond1F1
const pochhammer = HypergeometricFunctions.pochhammer

function w(z, par)
    @unpack w₀, z_R = par

    w₀ * √(1 + (z/z_R)^2)
end

R(z, z_R) = z + z_R^2 / z

include("envelopes.jl")
include("electricfield.jl")
include("magneticfield.jl")
include("potential.jl")
include("gauss.jl")
include("laguerre-gauss.jl")
include("setup.jl")

end # module
