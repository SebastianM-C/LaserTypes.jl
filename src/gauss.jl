# # Gaussian profile

struct GaussLaserConstantCache{IC,W,K,T,Z,E}
    inv_c::IC
    ω::W
    k::K
    T₀::T
    z_R::Z
    E₀::E
end

function GaussLaserConstantCache(;c, λ, w₀, a₀, mₑ, q)
    ω = 2π * c / λ
    k = 2π / λ
    T = 2π / ω
    z_R = w₀^2 * k / 2
    E₀ = a₀ * mₑ * c * ω / abs(q)

    return GaussLaserConstantCache(inv(c), ω, k, T, z_R, E₀)
end

@auto_hash_equals mutable struct GaussLaserCache{X,CE}
    x::X
    y::X
    wz::X
    Ex::CE
    Ey::CE
    Ez::CE
end

function GaussLaserCache(λ, E)
    TaskLocalValue{GaussLaserCache{typeof(zero(λ)),typeof(zero(E*im))}}() do
        GaussLaserCache(
            zero(λ),        # x
            zero(λ),        # y
            zero(λ),        # wz
            zero(E*im),     # Ex
            zero(E*im),     # Ey
            zero(E*im)      # Ez
        )
    end
end

function Base.fill!(cache::GaussLaserCache, x::AbstractVector)
    cache.x = x[1]
    cache.y = x[2]

    return nothing
end

@doc """
    struct GaussLaser <: AbstractLaser

The `GaussLaser` is defined by the `units` of the laser
(a positional argument that can be `:SI`, `:atomic` or `:SI_unitful` and `:atomic_unitful`)
and the following parameters
- `λ` is the laser wavelangth
- `a₀` is the normalized vector potential (defined as ``a_0=\\frac{eA}{m_e c^2}``)
- `ϕ₀` is the initial phase with the default value 0.0
- `w₀` is the beam radius at the Rayleigh range or [beam waist](https://en.wikipedia.org/wiki/Gaussian_beam#Beam_waist)
- `ξx` and `ξy` give the polarization and have the default value `1.0 + 0im` and `0.0 + 0im`
- `orientation` specifies how the laser is oriented with respect to the
default coordinate system (default `(:x, :z)`). See [`LaserGeometry`](@ref) for more details.
- `propagation_dir` is the propagation direction of the wave in the
intrinsic coordinate system (default `:z`). See [`LaserGeometry`](@ref) for more details.
- `profile` is the temporal profile of the pulse and the default one is constant (infinite pules duration)

During the initialization of the `GaussLaser` type, some useful derived values are
also computed
- `ω` is the angular frequency and is given by ``2π c / λ``
- `k` is the wavenumber and is given by ``2π / λ``
- `z_R` is the [Rayleigh range](https://en.wikipedia.org/wiki/Rayleigh_length) and is given by ``k w_0^2 / 2``
- `T₀` is the laser period and is given by ``2π / ω``
- `E₀` is the amplitude of the electric field and is given by ``a_0\\frac{m_q c \\omega}{q}``
"""
GaussLaser

struct GaussLaser{PhysicsConst,
                  Derived,
                  Cache,
                  Geometry,
                  Polarization,
                  Profile,
                  L,
                  F} <: AbstractLaser
    constants::PhysicsConst
    derived::Derived
    cache::Cache
    geometry::Geometry
    polarization::Polarization
    profile::Profile
    # laser parameters
    λ::L
    a₀::F
    ϕ₀::F
    w₀::L
end

function GaussLaser(units;
        λ,
        a₀,
        ϕ₀ = 0.0,
        w₀,
        ξx = 1.0+0im,
        ξy = 0,
        orientation = (:x, :z),
        profile = ConstantProfile()
    )

    ξx, ξy = promote(ξx, ξy)
    a₀, ϕ₀ = promote(a₀, ϕ₀)
    λ, w₀ = promote(λ, w₀)
    @assert hypot(ξx, ξy) ≈ 1 "Invalid ξx and ξy"
    constants = FundamentalConstants(units)
    @unpack mₑ, c, q = constants

    derived = GaussLaserConstantCache(; c, λ, w₀, a₀, mₑ, q)
    E₀ = derived.E₀

    cache = GaussLaserCache(λ, E₀)

    geometry = LaserGeometry(orientation)

    polarization = LaserPolarization(ξx, ξy)

    return GaussLaser(
        constants,
        derived,
        cache,
        geometry,
        polarization,
        profile,
        λ,
        a₀,
        ϕ₀,
        w₀
    )
end

update_cache!(::GaussLaser{P,D,Cache}, ::Any, ::Any) where {P,D,Cache<:Nothing} = nothing

function required_coords(::GaussLaser, r)
    CylindricalFromCartesian()(r)
end

function Ex(laser::GaussLaser, coords)
    @unpack E₀, k, z_R = immutable_cache(laser)
    ξx = polarization(laser, :ξx)
    @unpack w₀, ϕ₀ = laser
    @unpack r, z = coords

    wz = w(z, laser)
    update_cache!(laser, :wz, wz)

    ξx * E₀ * w₀/wz * exp(-(r/wz)^2 + im*(-(r^2*z)/(z_R*wz^2) + atan(z, z_R) - k*z + ϕ₀))
end

function Ez(laser::GaussLaser, coords)
    @unpack k, z_R = immutable_cache(laser)
    @unpack Ex, Ey, x, y, wz = mutable_cache(laser)
    z = coords.z

    2im/(k*wz^2)*(1 + im*(z/z_R))*(x*Ex + y*Ey)
end

function Bz(laser::GaussLaser, coords)
    @unpack k, z_R = immutable_cache(laser)
    @unpack Ex, Ey, x, y, wz = mutable_cache(laser)
    c = fundamental_constants(laser, :c)
    z = coords.z

    2im/(k*c*wz^2)*(1 + im*(z/z_R))*(y*Ex - x*Ey)
end
