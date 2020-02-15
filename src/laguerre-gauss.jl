using GLSL
using Unitful
using Parameters
using GeometryTypes: Vec3
import PhysicalConstants.CODATA2018: c_0, m_e, e
using ..LaserTypes: TemporalProfiles, w

@with_kw struct LaserParams{V,Q,M,L,T,F,C,I,P,W,K,E}
    # independent values
    c::V = c_0
    q::Q = -e
    m_q::M = m_e
    λ₀::L = 0.8u"μm"
    a₀::F = 1.0
    φ₀::F = 0.0
    w₀::L = 58.0u"μm"
    ξx::C = 1.0 + 0im
    ξy::C = -1.0im
    @assert hypot(ξx, ξy) ≈ 1
    l::I = 0
    p::I = 1
    τ₀::T = 18.02u"fs"
    z_F::L = uconvert(unit(λ₀), -4*τ₀*c)
    envelope::P = TemporalProfiles.gaussian
    # derived values
    ω₀::W = 2π * c / λ₀; @assert ω₀ ≈ 2π * c / λ₀
    k₀::K = 2π / λ₀; @assert k₀ ≈ 2π / λ₀
    z_R::L = w₀^2 * k₀ / 2; @assert z_R ≈ w₀^2 * k₀ / 2
    T₀::T = uconvert(unit(τ₀), 2π / ω₀); @assert T₀ ≈ 2π / ω₀
    E₀::E = a₀ * m_q * c * ω₀ / abs(q); @assert E₀ ≈ a₀ * m_q * c * ω₀ / abs(q)
end
