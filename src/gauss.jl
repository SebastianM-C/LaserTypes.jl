# # Gaussian profile

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

@with_kw struct GaussLaser{V,Q,M,L,F,C,T,P,W,K,E}
    # independent values
    c::V = c_0
    q::Q = -e
    m_q::M = m_e
    λ::L = 0.8u"μm"
    a₀::F = 1.0
    ϕ₀::F = 0.0
    w₀::L = 58.0u"μm"
    ξx::C = 1.0 + 0im
    ξy::C = 0.0 + 0im
    @assert hypot(ξx, ξy) ≈ 1
    profile::P = GaussianProfile(c=c)
    # derived values
    ω::W = 2π * c / λ; @assert ω ≈ 2π * c / λ
    k::K = 2π / λ; @assert k ≈ 2π / λ
    z_R::L = w₀^2 * k / 2; @assert z_R ≈ w₀^2 * k / 2
    T₀::T = 2π / ω; @assert T₀ ≈ 2π / ω
    E₀::E = a₀ * m_q * c * ω / abs(q); @assert E₀ ≈ a₀ * m_q * c * ω / abs(q)
end

function required_coords(laser::GaussLaser, r)
    PolarFromCartesian()(Vec2(r[1], r[2]))
end

function Ex(laser::GaussLaser, coords)
    @unpack E₀, w₀, k, z_R, ϕ₀, ξx = laser
    @unpack r, z = coords

    wz = w(z, laser)
    Rz = R(z, z_R)

    ξx * E₀ * w₀/wz * exp(-im*k*z - (r/wz)^2 - im*((k*r^2)/(2Rz) - atan(z, z_R) - ϕ₀))
end

function Ez(laser::GaussLaser, coords, Ex, Ey, x, y)
    @unpack k, z_R = laser
    z = coords.z

    wz = w(z, laser)

    2(im - z/z_R) / (k*wz^2) * (x*Ex + y*Ey)
end

function Bz(laser::GaussLaser, coords, Ex, Ey, x, y)
    @unpack k, z_R, c = laser
    z = coords.z

    wz = w(z, laser)

    2(im - z/z_R) / (c*k*wz^2) * (y*Ex - x*Ey)
end
