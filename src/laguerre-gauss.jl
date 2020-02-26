module LaGuerreGauss

using HypergeometricFunctions
using Unitful
using Parameters
using GeometryTypes: Vec3
import PhysicalConstants.CODATA2018: c_0, m_e, e
using ..LaserTypes: TemporalProfiles, w, g
using ..Gauss

const ₁F₁ = drummond1F1

@with_kw struct LaserParams{V,Q,M,L,F,C,T,P,I,W,K,E,R}
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
    Nₚₘ::R = √(pochhammer(p+1, abs(l))); @assert Nₚₘ ≈ √(pochhammer(p+1, abs(l)))
end

function Ex(x, y, z, r, par)
    @unpack Nₚₘ, w₀, z_R, ξx, p, m = par
    wz = w(z, par)
    Rz = R(z, z_R)
    Eg = Gauss.Ex(z, r, par)
    σ = (r/wz)^2
    mₐ = abs(m)
    φ = atan(x, y)

    ξx*Eg*Nₚₘ*(r*√2/wz)^mₐ*₁F₁(-p, mₐ+1, 2σ)*exp(im*((2p+mₐ)*atan(z, z_R)-m*φ-φ₀))
end

Ey(x, y, z, r, par) = par.ξy / par.ξx * Ex(x, y, z, r, par)

function Ez(x, y, z, r, par)
    @unpack k, z_R = par
    wz = w(z, par)

    -im / k * (-2*(1+im*(z/z_R))/wz + ∂₁F₁())
end

end  # module LaGuerreGauss
