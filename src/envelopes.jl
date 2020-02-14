module TemporalProfiles

using Parameters

function gaussian(z, t, par)
    @unpack ω₀, z_F, c, τ₀ = par

    exp(im*ω₀*t - ((t - (z-z_F) / c) / τ₀)^2)
end

end  # module TemporalProfiles
