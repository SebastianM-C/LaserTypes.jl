module LaGuerreGauss

using GSL: sf_laguerre_n
using FastTransforms: pochhammer
using Unitful
using Parameters
using GeometryTypes: Vec3
import PhysicalConstants.CODATA2018: c_0, m_e, e
using ..LaserTypes: TemporalProfiles, w, g

@with_kw struct LaserParams{V,Q,M,L,F,C,T,P,I,W,K,E,R}
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
    l::I = 0
    # derived values
    ω::W = 2π * c / λ; @assert ω ≈ 2π * c / λ
    k::K = 2π / λ; @assert k ≈ 2π / λ
    z_R::L = w₀^2 * k / 2; @assert z_R ≈ w₀^2 * k / 2
    T₀::T = uconvert(unit(τ₀), 2π / ω); @assert T₀ ≈ 2π / ω
    E₀::E = a₀ * m_q * c * ω / abs(q); @assert E₀ ≈ a₀ * m_q * c * ω / abs(q)
    Cₚₗ::R = 1 / √(pochhammer(p+1, abs(l)))
end

end  # module LaGuerreGauss
