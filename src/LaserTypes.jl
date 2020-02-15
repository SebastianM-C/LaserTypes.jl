module LaserTypes

export Gauss, TemporalProfiles

using Parameters

function w(z, par)
    @unpack w₀, z_R = par

    w₀ * √(1 + (z/z_R)^2)
end

R(z, z_R) = z + z_R^2 / z

"""
    g(z, t, par)

The time dependent part of the laser pulse.
"""
function g(z, t, par)
    @unpack envelope, ω = par

    exp(im*ω*t) * envelope(z, t, par)
end

include("envelopes.jl")
include("gauss.jl")
include("laguerre-gauss.jl")

# https://github.com/JuliaLang/julia/issues/34771
Base.hypot(x::T, y::T) where {T<:Number} = !iszero(x) ? (z = y/x; abs(x) * sqrt(one(z) + z*z)) : abs(y)

end # module
