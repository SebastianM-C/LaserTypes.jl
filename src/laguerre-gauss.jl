# # Laguerre-Gauss profile

struct LaguerreGaussLaserConstantCache{IC,W,K,T,Z,E,I,N}
    inv_c::IC
    ω::W
    k::K
    T₀::T
    z_R::Z
    E₀::E
    mₐ::I
    Nₚₘ::N
end

function LaguerreGaussLaserConstantCache(;c, λ, w₀, a₀, mₑ, q, p, m)
    ω = 2π * c / λ
    k = 2π / λ
    T = 2π / ω
    z_R = w₀^2 * k / 2
    E₀ = a₀ * mₑ * c * ω / abs(q)
    mₐ = abs(m)
    Nₚₘ = √(pochhammer(p+1, mₐ))

    return LaguerreGaussLaserConstantCache(inv(c), ω, k, T, z_R, E₀, mₐ, Nₚₘ)
end

@auto_hash_equals mutable struct LaguerreGaussLaserCache{L,S,CE,EE}
    x::L
    y::L
    σ::S
    wz::L
    Ex::CE
    Ey::CE
    Ez::CE
    Eg::CE
    NEgexp::EE
    rwz::S
end

function LaguerreGaussLaserCache(λ, E)
    ThreadLocal(LaguerreGaussLaserCache(
        zero(λ),                    # x
        zero(λ),                    # y
        zero(λ/λ),                  # σ
        zero(λ),                    # wz
        zero(E*im),                 # Ex
        zero(E*im),                 # Ey
        zero(E*im),                 # Ez
        zero(E*im),                 # Eg
        zero(E*im),                 # NEgexp
        zero(λ/λ),                  # rwz
    ))
end


function Base.fill!(cache::LaguerreGaussLaserCache, x::AbstractVector)
    cache.x = x[1]
    cache.y = x[2]

    return nothing
end

@doc """
    struct LaguerreGaussLaser <: AbstractLaser

The `LaguerreGaussLaser` is defined by the `units` of the laser
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
- `p` is the radial index of the mode, ``p ∈ ℤ, p ≥ 0``, with the default value 1
- `m` is the azimuthal index of the mode, ``m ∈ ℤ``, with the default value 0

During the initialization of the `LaguerreGaussLaser` type, some useful derived values are
also computed
- `ω` is the angular frequency and is given by ``2π c / λ``
- `k` is the wavenumber and is given by ``2π / λ``
- `z_R` is the [Rayleigh range](https://en.wikipedia.org/wiki/Rayleigh_length) and is given by ``k w_0^2 / 2``
- `T₀` is the laser period and is given by ``2π / ω``
- `E₀` is the amplitude of the electric field and is given by ``a_0\\frac{m_q c \\omega}{q}``
- `Nₚₘ` is a normalization factor given by ``\\sqrt{(p+1)_{|m|}}``, with ``(x)_n`` the [Pochhammer symbol](https://mathworld.wolfram.com/PochhammerSymbol.html)
"""
LaguerreGaussLaser

struct LaguerreGaussLaser{_P,_M,
                          C0,Q,M,Eps,Mu,U,
                          IC,W,K,T,Z,E,I,F,
                          L,S,CE,EE,
                          D,R,
                          C,
                          P} <: AbstractLaser
    constants::FundamentalConstants{C0,Q,M,Eps,Mu,U}
    derived::LaguerreGaussLaserConstantCache{IC,W,K,T,Z,E,I,F}
    cache::ThreadLocal{LaguerreGaussLaserCache{L,S,CE,EE}}
    geometry::LaserGeometry{D,R}
    polarization::LaserPolarization{C}
    profile::P
    # laser parameters
    λ::L
    a₀::F
    ϕ₀::F
    w₀::L
    p::I
    m::I
end

# This hack stores the m and p values in the type domain
function LaguerreGaussLaser(
    constants::FundamentalConstants{C0,Q,M,Eps,Mu,U},
    derived::LaguerreGaussLaserConstantCache{IC,W,K,T,Z,E,I,F},
    cache::ThreadLocal{LaguerreGaussLaserCache{L,S,CE,EE}},
    geometry::LaserGeometry{D,R},
    polarization::LaserPolarization{C},
    profile::P,
    λ::L,
    a₀::F,
    ϕ₀::F,
    w₀::L,
    p::I,
    m::I
) where {C0,Q,M,Eps,Mu,U,
         IC,W,K,T,Z,E,I,F,
         L,S,CE,EE,
         D,R,
         C,
         P}
    LaguerreGaussLaser{Val{p},Val{m},
                       C0,Q,M,Eps,Mu,U,
                       IC,W,K,T,Z,E,I,F,
                       L,S,CE,EE,
                       D,R,
                       C,
                       P}(
        constants,
        derived,
        cache,
        geometry,
        polarization,
        profile,
        λ,
        a₀,
        ϕ₀,
        w₀,
        p,
        m
    )
end

function LaguerreGaussLaser(units;
        λ,
        a₀,
        ϕ₀ = 0.0,
        w₀,
        p = 1,
        m = 0,
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

    derived = LaguerreGaussLaserConstantCache(; c, λ, w₀, a₀, mₑ, q, p, m)
    E₀ = derived.E₀

    cache = LaguerreGaussLaserCache(λ, E₀)

    geometry = LaserGeometry(orientation)

    polarization = LaserPolarization(ξx, ξy)

    return LaguerreGaussLaser(
        constants,
        derived,
        cache,
        geometry,
        polarization,
        profile,
        λ,
        a₀,
        ϕ₀,
        w₀,
        p,
        m
    )
end

"""
    convert_laser(::Type{GaussLaser}, laser::LaguerreGaussLaser)

Convert a `LaguerreGaussLaser` to a `GaussLaser` with similar parameters.
"""
function convert_laser(::Type{GaussLaser}, laser::LaguerreGaussLaser)
    @unpack constants, geometry, profile, λ, a₀, w₀ = laser

    derived = GaussLaserConstantCache(;
        constants.c,
        λ,
        w₀,
        a₀,
        constants.mₑ,
        constants.q)
    # E₀ = derived.E₀

    # Don't create a cache since it's not used and it allocates memory
    cache = nothing # GaussLaserCache(λ, E₀)

    polarization = LaserPolarization(1, 0)

    GaussLaser(
        constants,
        derived,
        cache,
        geometry,
        polarization,
        profile,
        λ,
        a₀,
        zero(a₀),
        w₀
    )
end

get_p(::Type{<:LaguerreGaussLaser{Val{P},Val{M}}}) where {P,M} = P
get_m(::Type{<:LaguerreGaussLaser{Val{P},Val{M}}}) where {P,M} = M

function required_coords(::LaguerreGaussLaser, r)
    CylindricalFromCartesian()(r)
end

function Ex(laser::LaguerreGaussLaser, coords)
    @unpack Nₚₘ, z_R, mₐ = immutable_cache(laser)
    ξx = polarization(laser, :ξx)
    @unpack ϕ₀, p, m, cache = laser
    @unpack r, θ, z = coords

    gauss_laser = convert_laser(GaussLaser, laser)
    Eg = Ex(gauss_laser, coords)
    # We avoid computing wz in the GaussLaser call so that we can skip cache creation
    wz = w(z, laser)
    σ = (r/wz)^2
    rwz = ustrip(NoUnits, r*√2/wz)
    @pack! cache[] = Eg, wz, σ, rwz

    ξx * Eg * Nₚₘ * rwz^mₐ * _₁F₁(-p, mₐ+1, 2σ) * exp(im*((2p+mₐ)*atan(z, z_R) - m*θ + ϕ₀))
end

function Ez(laser::LaguerreGaussLaser, coords)
    @unpack Nₚₘ, mₐ, k, z_R = immutable_cache(laser)
    @unpack wz, rwz, σ, Eg, Ex, Ey, x, y = mutable_cache(laser)
    @unpack ξx, ξy = polarization(laser)
    @unpack ϕ₀, p, m = laser
    @unpack r, θ, z = coords

    P = get_p(typeof(laser))
    M = get_m(typeof(laser))

    sgn = sign(m)
    𝟘 = zero(typeof(Ex))/oneunit(typeof(x))
    NEgexp = Nₚₘ*Eg*exp(im*((2p+mₐ)*atan(z, z_R)-m*θ+ϕ₀))
    update_cache!(laser, :NEgexp, NEgexp)

    -im/k * (
       (iszero(M) ? 𝟘 : mₐ * (ξx - im*sgn*ξy) * (√2/wz)^mₐ * r^(mₐ-1) * _₁F₁(-p, mₐ+1, 2σ) * NEgexp * exp(im*sgn*θ))
     - 2/(wz^2) * (1 + im*z/z_R) * (x*Ex + y*Ey)
     - (iszero(P) ? 𝟘 : 4p/((mₐ+1) * wz^2) * (x*ξx + y*ξy) * rwz^mₐ * _₁F₁(-p+1, mₐ+2, 2σ) * NEgexp)
    )
end

function Bz(laser::LaguerreGaussLaser, coords)
    @unpack ω, z_R = immutable_cache(laser)
    @unpack wz, mₐ, σ, Ex, Ey, NEgexp, x, y = mutable_cache(laser)
    @unpack ξx, ξy = polarization(laser)
    @unpack p, m = laser
    @unpack r, θ, z = coords

    sgn = sign(m)
    𝟘 = zero(typeof(Ex))/oneunit(typeof(x))

    -im/ω * (
       - (iszero(m) ? 𝟘 : mₐ * (ξy + im*sgn*ξx) * (√2/wz)^mₐ*r^(mₐ-1) * _₁F₁(-p, mₐ+1, 2σ) * NEgexp * exp(im*sgn*θ))
       + 2/(wz^2) * (1 + im*z/z_R) * (x*Ey - y*Ex)
       + (iszero(p) ? 𝟘 : (4p)/((mₐ+1) * wz^2) * (x*ξy - y*ξx) * (r*√2/wz)^mₐ * _₁F₁(-p+1, mₐ+2, 2σ) * NEgexp)
    )
end
