# # Envelopes

"""
    struct ConstantProfile

This is the trivial profile
```math
envelope(z, t) = 1
```
which gives an infinite duration for the electromagnetic field configuration (we cannot
call this a laser pusle since a pulse implicitely has a finite duration).
"""
struct ConstantProfile end

@doc """
    GaussianProfile{V,T,L}

This envelope provides a finite duration for the laser pulse and thus can provide a more
realistic description of an actual laser pulse.
```math
envelope(z, t) = \\cosh\\left[\\left(\\frac{\\varphi}{\\tau})\\right)^2\\right],
```
where
```math
\\varphi = (t - t_0) - \\frac{z - z_0}{c} ,
```
and
- `c` is the speed of light
- `τ` is the duration of the pulse (FWHM) and has the default value 18.02fs
- `t₀` is the origin of the time axis and it is 0 by default
- `z₀` is the initial position of the intensity peak and has the default value `-4*τ*c`
"""
GaussianProfile

@with_kw struct GaussianProfile{V,T,L}
    c::V = c_0
    τ::T = 18.02u"fs"
    t₀::T = zero(τ)
    z₀::L = -4*τ*c
end

@doc """
    QuasiRectangularProfile{V,T,L}

This envelope provides a finite duration for the laser pulse and thus can provide a more
realistic description of an actual laser pulse.
```math
envelope(z, t) =
```
where
- `c` is the speed of light
- `τ` is the duration of the pulse (FWHM) and has the default value 18.02fs
- `t₀` is the origin of the time axis and it is 0 by default
- `z₀` is the initial position of the intensity peak and has the default value `-4*τ*c`
- `Δz` is the width of the flat part of the profile and the default value `10*τ*c`
"""
QuasiRectangularProfile

@with_kw struct QuasiRectangularProfile{V,T,L}
    c::V = c_0
    τ::T = 18.02u"fs"
    t₀::T = zero(τ)
    z₀::L = -4*τ*c
    Δz::L = 10*τ*c
end

"""
    g(z, t, par)

The time dependence of the fields is given by
```math
g(z, t) = \\cos(\\omega t) envelope(z, t),
```
where
- ``\\omega`` is the angular frequency of the laser pulse
- ``envelope(z, t)`` is a function that can be used to control the duration of the pulse

"""
function g(z, t, par)
    @unpack profile, ω = par

    cos(ω*t) * envelope(profile, z, t)
end

envelope(::ConstantProfile, z, t) = 1

function envelope(profile::GaussianProfile, z, t)
    @unpack c, τ, t₀, z₀ = profile
    φ = (t - t₀) - (z - z₀) / c

    cosh(uconvert(NoUnits, φ / τ)^2)
end

function envelope(profile::QuasiRectangularProfile, z, t)
    @unpack c, τ, t₀, z₀, Δz = profile
    φ = (t - t₀) - (z - z₀) / c

    if φ < zero(φ)
        cosh(uconvert(NoUnits, φ / τ)^2)
    elseif φ < Δz
        one(φ)
    else
        cosh(uconvert(NoUnits, (φ - Δz/c) / τ)^2)
    end
end
