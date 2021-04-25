# # Laguerre-Gauss profile

@doc """
    struct LaguerreGaussLaser{V,Q,M,L,F,C,T,P,I,W,K,E,R}

The `LaguerreGaussLaser` is defined by the following independent parameters
- `c` is the speed of light in vaccum, with the default value being in SI (`c_0` from the CODATA2018 in the [PhysicalConstants](https://github.com/JuliaPhysics/PhysicalConstants.jl) package)
- `q` is the electric charge, with the default value being the one for the electron in SI (`-e` from the CODATA2018 in the [PhysicalConstants](https://github.com/JuliaPhysics/PhysicalConstants.jl) package)
- `m_q` is the mass of the charge, with the default value being the one for the electron in SI (`m_e` from the CODATA2018 in the [PhysicalConstants](https://github.com/JuliaPhysics/PhysicalConstants.jl) package)
- `λ` is the laser wavelangth with the default value 0.8μm
- `a₀` is the normalized vector potential (defined as ``a_0=\\frac{eA}{m_e c^2}``)
- `ϕ₀` is the initial phase with the default value 0.0
- `w₀` is the beam radius at the Rayleigh range or [beam waist](https://en.wikipedia.org/wiki/Gaussian_beam#Beam_waist) with the default value 58.0μm
- `ξx` and `ξy` give the polarization and have the default value `1.0 + 0im` and `0.0 + 0im`
- `profile` is the temporal profile of the pulse and the default one is a Gaussian one
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

@with_kw struct LaguerreGaussLaser{V,Q,M,M0,L1,L2,L3,F,C,T,P,I,W,K,E,R}
    # independent values
    c::V = c_0
    q::Q = -e
    m_q::M = m_e
    μ₀::M0 = μ_0
    λ::L1 = 0.8u"μm"
    a₀::F = 1.0
    ϕ₀::F = 0.0
    w₀::L2 = 58.0u"μm"
    ξx::C = 1.0 + 0im
    ξy::C = 0.0 + 0im
    @assert hypot(ξx, ξy) ≈ 1
    profile::P = GaussProfile(c=c)
    p::I = 1
    m::I = 0
    # derived values
    ω::W = 2π * c / λ; @assert ω ≈ 2π * c / λ
    k::K = 2π / λ; @assert k ≈ 2π / λ
    z_R::L3 = w₀^2 * k / 2; @assert z_R ≈ w₀^2 * k / 2
    T₀::T = 2π / ω; @assert T₀ ≈ 2π / ω
    E₀::E = a₀ * m_q * c * ω / abs(q); @assert E₀ ≈ a₀ * m_q * c * ω / abs(q)
    Nₚₘ::R = √(pochhammer(p+1, abs(m))); @assert Nₚₘ ≈ √(pochhammer(p+1, abs(m)))
end

"""
    Base.convert(::Type{GaussLaser}, laser::LaguerreGaussLaser)

Convert a `LaguerreGaussLaser` to a `GaussLaser` with the same parameters.
"""
function Base.convert(::Type{GaussLaser}, laser::LaguerreGaussLaser)
    @unpack c, q, m_q, μ₀, λ, a₀, ϕ₀, w₀, ξx, ξy, profile, ω, k, z_R, T₀, E₀ = laser
    GaussLaser(c, q, m_q, μ₀, λ, a₀, ϕ₀, w₀, ξx, ξy, profile, ω, k, z_R, T₀, E₀)
end

function required_coords(laser::LaguerreGaussLaser, r)
    CylindricalFromCartesian()(r)
end

function Ex(laser::LaguerreGaussLaser, coords)
    @unpack Nₚₘ, w₀, ϕ₀, z_R, ξx, p, m = laser
    @unpack r, θ, z = coords

    wz = w(z, laser)
    gauss_laser = convert(GaussLaser, laser)
    Eg = Ex(gauss_laser, coords)
    σ = (r/wz)^2
    mₐ = abs(m)

    ξx*Eg*Nₚₘ*(r*√2/wz)^mₐ*_₁F₁(-p, mₐ+1, 2σ)*exp(im*((2p+mₐ)*atan(z, z_R)-m*θ-ϕ₀))
end

function Ez(laser::LaguerreGaussLaser, coords, E_x, E_y, x, y)
    @unpack Nₚₘ, w₀, ϕ₀, k, z_R, p, m, ξx, ξy = laser
    @unpack r, θ, z = coords

    wz = w(z, laser)
    mₐ = abs(m)
    σ = (r/wz)^2
    ∓ = m > 0 ? (-) : +
    gauss_laser = convert(GaussLaser, laser)
    Eg = Ex(gauss_laser, coords)
    𝟘 = zero(typeof(E_x))/oneunit(typeof(x))
    ExImEy = Eg*Nₚₘ*r^(mₐ-1)*(√2/wz)^mₐ*_₁F₁(-p, mₐ+1, 2σ)*exp(im*((2p+mₐ)*atan(z, z_R)-(m+1)*θ-ϕ₀))

    -im / k * (
        -2*(1+im*(z/z_R))/wz^2 * (x*E_x + y*E_y)
        + 4p/((mₐ+1)*wz^2) * (x*ξx+y*ξy) * Eg*Nₚₘ*(r*√2/wz)^mₐ*exp(im*((2p+mₐ)*atan(z, z_R)-m*θ-ϕ₀))  
        - (!iszero(m) ? mₐ * (ξx ∓ im*ξy) * ExImEy : 𝟘)
        )
end

function Bz(laser::LaguerreGaussLaser, coords, E_x, E_y, x, y)
    @unpack Nₚₘ, w₀, ϕ₀, k, c, z_R, p, m, ξx, ξy = laser
    @unpack r, θ, z = coords

    wz = w(z, laser)
    σ = (r/wz)^2
    mₐ = abs(m)
    ∓ = m > 0 ? (-) : +
    gauss_laser = convert(GaussLaser, laser)
    Eg = Ex(gauss_laser, coords)
    𝟘 = zero(typeof(E_x))/oneunit(typeof(x))
    ExImEy = Eg*Nₚₘ*r^(mₐ-1)*(√2/wz)^mₐ*_₁F₁(-p, mₐ+1, 2σ)*exp(im*((2p+mₐ)*atan(z, z_R)-(m+1)*θ-ϕ₀))

    -im / (k*c) * (
        -2*(1+im*(z/z_R))/wz^2 * (x*E_x + y*E_y)
        + 4p/((mₐ+1)*wz^2) * (x*ξy+y*ξx) * Eg*Nₚₘ*(r*√2/wz)^mₐ*exp(im*((2p+mₐ)*atan(z, z_R)-m*θ-ϕ₀))  
        - (!iszero(m) ? mₐ * (ξx ∓ im*ξy) * ExImEy : 𝟘)
        )
end
