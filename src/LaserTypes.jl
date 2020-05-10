# # LaserTypes

module LaserTypes

export E, B, GaussLaser, LaguerreGaussLaser, TemporalProfiles

using Parameters

using Unitful
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

"""
    g(z, t, par)

The time dependence of the fields is given by
```math
g(z, t) = \\exp(\\mathrm{i} \\omega t) envelope(z, t),
```
where
- ``\\omega`` is the angular frequency of the laser pulse
- ``envelope(z, t)`` is a function that can be used to control the duration of the pulse

"""
function g(z, t, par)
    @unpack envelope, ω = par

    exp(im*ω*t) * envelope(z, t, par)
end

include("envelopes.jl")
include("electricfield.jl")
include("magneticfield.jl")
include("potential.jl")
include("gauss.jl")
include("laguerre-gauss.jl")

end # module
