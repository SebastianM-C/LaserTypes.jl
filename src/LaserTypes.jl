module LaserTypes

export Gauss, TemporalProfiles

using Parameters

function w(z, par)
    @unpack w₀, z_R = par

    w₀ * √(1 + (z/z_R)^2)
end

function g(z, t, par)
    @unpack envelope, ω = par

    exp(im*ω*t) * envelope(z, t, par)
end

include("envelopes.jl")
include("gauss.jl")


end # module
