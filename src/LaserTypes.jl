# # LaserTypes

module LaserTypes

export E, B, GaussLaser, LaguerreGaussLaser, ConstantProfile, GaussianProfile, QuasiRectangularProfile

using Unitful
using Parameters
using HypergeometricFunctions
using GeometryTypes: Vec2, Vec3
using CoordinateTransformations
import PhysicalConstants.CODATA2018: c_0, m_e, e

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

end # module
