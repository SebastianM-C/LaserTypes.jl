# # Laguerre-Gauss profile

@doc """
    struct LaguerreGaussLaser{V,Q,M,L,F,C,T,P,I,W,K,E,R}

The `LaguerreGaussLaser` is defined by the following independent parameters
- `c` is the speed of light in vaccum, with the default value being in SI (`c_0` from the CODATA2018 in the [PhysicalConstants](https://github.com/JuliaPhysics/PhysicalConstants.jl) package)
- `q` is the electric charge, with the default value being the one for the electron in SI (`-e` from the CODATA2018 in the [PhysicalConstants](https://github.com/JuliaPhysics/PhysicalConstants.jl) package)
- `m_q` is the mass of the charge, with the default value being the one for the electron in SI (`m_e` from the CODATA2018 in the [PhysicalConstants](https://github.com/JuliaPhysics/PhysicalConstants.jl) package)
- `λ` is the laser wavelangth with the default value 0.8μm
- `a₀` is the normalized vector potential (defined as ``a_0=\\frac{eA}{m_e c^2}``)
- `φ₀` is the initial phase with the default value 0.0
- `w₀` is the beam radius at the Rayleigh range or [beam waist](https://en.wikipedia.org/wiki/Gaussian_beam#Beam_waist) with the default value 58.0μm
- `ξx` and `ξy` give the polarization and have the default value `1.0 + 0im` and `0.0 + 0im`
- `envelope` is the temporal envelope of the pulse and the default one is a Gaussian one
- `p` with the default value 1
- `m` with the default value 0

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

@with_kw struct LaguerreGaussLaser{V,Q,M,L,F,C,T,P,I,W,K,E,R}
    # independent values
    c::V = c_0
    q::Q = -e
    m_q::M = m_e
    λ::L = 0.8u"μm"
    a₀::F = 1.0
    φ₀::F = 0.0
    w₀::L = 58.0u"μm"
    ξx::C = 1.0 + 0im
    ξy::C = 0.0 + 0im
    @assert hypot(ξx, ξy) ≈ 1
    τ₀::T = 18.02u"fs"
    z_F::L = uconvert(unit(λ), -4*τ₀*c)
    envelope::P = TemporalProfiles.gaussian
    p::I = 1
    m::I = 0
    # derived values
    ω::W = 2π * c / λ; @assert ω ≈ 2π * c / λ
    k::K = 2π / λ; @assert k ≈ 2π / λ
    z_R::L = w₀^2 * k / 2; @assert z_R ≈ w₀^2 * k / 2
    T₀::T = uconvert(unit(τ₀), 2π / ω); @assert T₀ ≈ 2π / ω
    E₀::E = a₀ * m_q * c * ω / abs(q); @assert E₀ ≈ a₀ * m_q * c * ω / abs(q)
    Nₚₘ::R = √(pochhammer(p+1, abs(m))); @assert Nₚₘ ≈ √(pochhammer(p+1, abs(m)))
end

"""
    Base.convert(::Type{GaussLaser}, laser::LaguerreGaussLaser)

Convert a `LaguerreGaussLaser` to a `GaussLaser` with the same parameters.
"""
function Base.convert(::Type{GaussLaser}, laser::LaguerreGaussLaser)
    @unpack c, q, m_q, λ, a₀, φ₀, w₀, ξx, ξy, τ₀, z_F, envelope, ω, k, z_R, T₀, E₀ = laser
    GaussLaser(c, q, m_q, λ, a₀, φ₀, w₀, ξx, ξy, τ₀, z_F, envelope, ω, k, z_R, T₀, E₀)
end

function Ex(laser::LaguerreGaussLaser, x, y, z, r)
    @unpack Nₚₘ, w₀, φ₀, z_R, ξx, p, m = laser
    wz = w(z, laser)
    Rz = R(z, z_R)
    gauss_laser = convert(GaussLaser, laser)
    Eg = Ex(gauss_laser, z, r)
    σ = (r/wz)^2
    mₐ = abs(m)
    φ = atan(x, y)

    ξx*Eg*Nₚₘ*(r*√2/wz)^mₐ*_₁F₁(-p, mₐ+1, 2σ)*exp(im*((2p+mₐ)*atan(z, z_R)-m*φ-φ₀))
end

function Ez(laser::LaguerreGaussLaser, Ex, Ey, x, y, z, r)
    @unpack k, z_R, p, m = laser
    wz = w(z, laser)
    mₐ = abs(m)
    ∓ = m > 0 ? (-) : +

    -im / k * (
        (-2*(1+im*(z/z_R))/wz^2
        + 4p/(((mₐ+1)*wz^2) * _₁F₁(-p+1, mₐ+2, 2r^2/wz^2))) * (x*Ex + y*Ey)
        - mₐ/(x+im*y) * (Ex ∓ im*Ey)
        )
end

function Bz(laser::LaguerreGaussLaser, Ex, Ey, x, y, z, r)
    @unpack k, z_R, p, m, c = laser
    wz = w(z, laser)
    mₐ = abs(m)
    ∓ = m > 0 ? (-) : +

    -im / (c*k) * (
        (-2*(1+im*(z/z_R))/wz^2
        + 4p/(((mₐ+1)*wz^2) * _₁F₁(-p+1, mₐ+2, 2r^2/wz^2))) * (y*Ex + x*Ey)
        - mₐ/(x+im*y) * (Ey ∓ im*Ex)
        )
end
