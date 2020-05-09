# # Laguerre-Gauss profile

@with_kw struct LaguerreGaussLaser{V,Q,M,L,F,C,T,P,I,W,K,E,R}
    # independent values:
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
    GaussLaser(laser::LaguerreGaussLaser)

Convert a `LaguerreGaussLaser` to a `GaussLaser` with the same parameters.
"""
function GaussLaser(laser::LaguerreGaussLaser)
    @unpack c, q, m_q, λ, a₀, φ₀, w₀, ξx, ξy, τ₀, z_F, envelope, ω, k, z_R, T₀, E₀ = laser
    GaussLaser(c, q, m_q, λ, a₀, φ₀, w₀, ξx, ξy, τ₀, z_F, envelope, ω, k, z_R, T₀, E₀)
end

function Ex(laser::LaguerreGaussLaser, x, y, z, r)
    @unpack Nₚₘ, w₀, φ₀, z_R, ξx, p, m = laser
    wz = w(z, laser)
    Rz = R(z, z_R)
    gauss_laser = GaussLaser(laser)
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
