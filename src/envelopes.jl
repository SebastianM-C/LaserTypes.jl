module TemporalProfiles

using Parameters

constant(z, t, par) = 1

function gaussian(z, t, par)
    @unpack c, τ₀, z_F = par

    exp(-((t - (z-z_F) / c) / τ₀)^2)
end

function gaussian(ϕ, par)
    @unpack c, τ₀ = par

    exp(-(ϕ / (c*τ₀))^2)
end

quasi_rectangular(z, t, par) = quasi_rectangular(par.c*t - z, par)

function quasi_rectangular(ϕ, par)
    @unpack τ₀, ϕc = par

    if ϕ < 0
        exp(-(ϕ / τ₀)^2)
    elseif ϕ < ϕc
        1
    else
        exp(-((ϕ - ϕc) / τ₀)^2)
    end
end

end  # module TemporalProfiles
