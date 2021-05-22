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
    GaussLaserCache(
        zero(λ),        # x
        zero(λ),        # y
        zero(λ),        # wz
        zero(E*im),     # Ex
        zero(E*im),     # Ey
        zero(E*im)      # Ez
    )
end

function Base.fill!(cache::GaussLaserCache, x::AbstractVector)
    cache.x = x[1]
    cache.y = x[2]

    return nothing
end

@doc """
    struct GaussLaser{V,Q,M,L,F,C,T,P,W,K,E}

The `GaussLaser` is defined by the following independent parameters
- `c` is the speed of light in vaccum, with the default value being in SI (`c_0` from the CODATA2018 in the [PhysicalConstants](https://github.com/JuliaPhysics/PhysicalConstants.jl) package)
- `q` is the electric charge, with the default value being the one for the electron in SI (`-e` from the CODATA2018 in the [PhysicalConstants](https://github.com/JuliaPhysics/PhysicalConstants.jl) package)
- `m_q` is the mass of the charge, with the default value being the one for the electron in SI (`m_e` from the CODATA2018 in the [PhysicalConstants](https://github.com/JuliaPhysics/PhysicalConstants.jl) package)
- `λ` is the laser wavelangth with the default value 0.8μm
- `a₀` is the normalized vector potential (defined as ``a_0=\\frac{eA}{m_e c^2}``)
- `ϕ₀` is the initial phase with the default value 0.0
- `w₀` is the beam radius at the Rayleigh range or [beam waist](https://en.wikipedia.org/wiki/Gaussian_beam#Beam_waist) with the default value 58.0μm
- `ξx` and `ξy` give the polarization and have the default value `1.0 + 0im` and `0.0 + 0im`
- `profile` is the temporal profile of the pulse and the default one is a Gaussian one

During the initialization of the `GaussLaser` type, some useful derived values are
also computed
- `ω` is the angular frequency and is given by ``2π c / λ``
- `k` is the wavenumber and is given by ``2π / λ``
- `z_R` is the [Rayleigh range](https://en.wikipedia.org/wiki/Rayleigh_length) and is given by ``k w_0^2 / 2``
- `T₀` is the laser period and is given by ``2π / ω``
- `E₀` is the amplitude of the electric field and is given by ``a_0\\frac{m_q c \\omega}{q}``
"""
GaussLaser

struct GaussLaser{C0,Q,M,Eps,Mu,IC,W,K,T,Z,E,L,CE,D,R,C,P,F} <: AbstractLaser
    constants::FundamentalConstants{C0,Q,M,Eps,Mu}
    derived::GaussLaserConstantCache{IC,W,K,T,Z,E}
    cache::GaussLaserCache{L,CE}
    geometry::LaserGeometry{D,R}
    polarization::LaserPolarization{C}
    profile::P
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
        oscillation_dir = :x,
        propagation_dir = :z,
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

    geometry = LaserGeometry(oscillation_dir, propagation_dir)

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
    Rz = R(z, z_R)

    ξx * E₀ * w₀/wz * exp(-im*k*z - (r/wz)^2 - im*((k*r^2)/(2Rz) - atan(z, z_R) - ϕ₀))
end

function Ez(laser::GaussLaser, coords)
    @unpack k, z_R = immutable_cache(laser)
    @unpack Ex, Ey, x, y, wz = mutable_cache(laser)
    z = coords.z

    2(im - z/z_R) / (k*wz^2) * (x*Ex + y*Ey)
end

function Bz(laser::GaussLaser, coords)
    @unpack k, z_R = immutable_cache(laser)
    @unpack Ex, Ey, x, y, wz = mutable_cache(laser)
    c = fundamental_constants(laser, :c)
    z = coords.z

    2(im - z/z_R) / (c*k*wz^2) * (y*Ex - x*Ey)
end
