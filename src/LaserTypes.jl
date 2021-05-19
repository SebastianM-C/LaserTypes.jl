# # LaserTypes

module LaserTypes

export E, B,
    GaussLaser, LaguerreGaussLaser,
    ConstantProfile, GaussProfile, Cos²Profile, QuasiRectangularProfile,
    setup_laser, Fμν, S

using Unitful
using UnitfulAtomic
using Parameters
using LinearAlgebra
using HypergeometricFunctions
using StaticArrays
using CoordinateTransformations
import PhysicalConstants.CODATA2018: c_0, m_e, e, μ_0

const _₁F₁ = HypergeometricFunctions.drummond1F1
const pochhammer = HypergeometricFunctions.pochhammer

function w(z, par)
    @unpack w₀, z_R = par

    w₀ * √(1 + (z/z_R)^2)
end

# R(z, z_R) = z + z_R^2 / z

include("envelopes.jl")
include("electricfield.jl")
include("magneticfield.jl")
include("faraday.jl")
include("gauss.jl")
include("laguerre-gauss.jl")
include("setup.jl")
include("derived.jl")

end # module
