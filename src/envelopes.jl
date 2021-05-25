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
struct ConstantProfile <: AbstractTemporalProfile end

@doc """
    GaussProfile{V,T,L}

This envelope provides a finite duration for the laser pulse and thus can provide a more
realistic description of an actual laser pulse.
```math
envelope(z, t) = \\exp\\left[\\left(\\frac{φ}{τ}\\right)^2\\right],
```
where
```math
\\varphi = (t - t_0) - \\frac{z - z_0}{c}\\,,
```
and
- `c` is the speed of light
- `τ` is the duration of the pulse (FWHM) and has the default value 18.02fs
- `t₀` is the origin of the time axis and it is 0 by default
- `z₀` is the initial position of the intensity peak and has the default value `-4*τ*c`
"""
GaussProfile

struct GaussProfile{T,IT,L} <: AbstractTemporalProfile
    inv_τ::IT
    t₀::T
    z₀::L
end

function GaussProfile(;τ, t₀=zero(τ), z₀)
    GaussProfile(inv(τ), t₀, z₀)
end

@doc """
    Cos²Profile{V,T,L}

This envelope provides a finite duration for the laser pulse and thus can provide a more
realistic description of an actual laser pulse.
```math
envelope(z, t) =
    \\begin{cases}
    \\cos\\left[\\left(\\frac{φ}{τ}\\right)\\right]^2, & \\text{for } |t - t₀| < τ / 2\\\\
    0 \\,, & \\text{otherwise}
    \\end{cases}
```
where
```math
\\varphi = (t - t_0) - \\frac{z - z_0}{c}\\,,
```
and
- `c` is the speed of light
- `τ` is the duration of the pulse and has the default value 18.02fs
- `t₀` is the origin of the time axis and it is 0 by default
- `z₀` is the initial position of the intensity peak and is 0 by default
"""
Cos²Profile

struct Cos²Profile{T,IT,L} <: AbstractTemporalProfile
    τ::T
    inv_τ::IT
    t₀::T
    z₀::L
end

function Cos²Profile(;τ, t₀=zero(τ), z₀)
    Cos²Profile(τ, inv(τ), t₀, z₀)
end

@doc """
    QuasiRectangularProfile{V,T,L}

This envelope produces a pulse with a predominant constant part of width ``Δz``
which could offer better results than the Gaussian profile in the paraxial limit
(which is considered for the spatial profiles). The shape of the envelope is given by
```math
envelope(z, t) =
    \\begin{cases}
    \\exp\\left[\\left(\\frac{φ}{τ}\\right)^2\\right], & \\text{for } φ ≤ 0\\\\
    1\\,, & \\text{for } φ < Δz\\\\
    \\exp\\left[\\left(\\frac{φ - Δt/c}{τ}\\right)^2\\right], & \\text{otherwise}
    \\end{cases}
```
where
```math
\\varphi = (t - t_0) - \\frac{z - z_0}{c}\\,,
```
and
- `c` is the speed of light
- `τ` is the duration of the pulse (FWHM) and has the default value 18.02fs
- `t₀` is the origin of the time axis and it is 0 by default
- `z₀` is the initial position of the intensity peak and has the default value `-4*τ*c`
- `Δt` is the duration of the flat part of the profile and the default value `10*τ`
"""
QuasiRectangularProfile

struct QuasiRectangularProfile{T,IT,L} <: AbstractTemporalProfile
    inv_τ::IT
    t₀::T
    z₀::L
    Δt::T
end

function QuasiRectangularProfile(;τ, t₀=zero(τ), z₀, Δt=10τ)
    t₀, Δt = promote(t₀, Δt)
    QuasiRectangularProfile(inv(τ), t₀, z₀, Δt)
end

"""
    g(z, t, par)

The time dependence of the fields defined by this package is given by
```math
g(z, t) = \\exp(i ω t) envelope(z, t),
```
where
- `z` and `t` are the position on the ``Oz`` axis and the time
- `par` are the laser parameters which pass the corresponding profile to the envelope
and
- ``ω`` is the angular frequency of the laser pulse
- ``envelope(z, t)`` is a function that can be used to control the duration of the pulse
"""
function g(z, t, laser; inv_c)
    profile = laser.profile
    ω = immutable_cache(laser, :ω)

    exp(im*ω*t) * envelope(profile, z, t; inv_c)
end

@inline envelope(::ConstantProfile, z, t; inv_c) = 1

@inline function envelope(profile::GaussProfile, z, t; inv_c)
    @unpack inv_τ, t₀, z₀ = profile
    φ = (t - t₀) - (z - z₀) * inv_c

    exp(-(φ * inv_τ)^2)
end

@inline function envelope(profile::Cos²Profile, z, t; inv_c)
    @unpack inv_τ, τ, t₀, z₀ = profile
    φ = π*((t - t₀) - (z - z₀) * inv_c)

    if abs(φ) / π < τ / 2
        cos(-(φ * inv_τ))^2
    else
        zero(t * inv_τ)
    end
end

@inline function envelope(profile::QuasiRectangularProfile, z, t; inv_c)
    @unpack inv_τ, t₀, z₀, Δt = profile
    φ = (t - t₀) - (z - z₀) * inv_c

    if φ < zero(φ)
        exp((φ * inv_τ)^2)
    elseif φ < Δt
        one(φ)
    else
        exp(((φ - Δz*inv_c) * inv_τ)^2)
    end
end
