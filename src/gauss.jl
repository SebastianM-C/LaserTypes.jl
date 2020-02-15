module Gauss

using Unitful
using Parameters
using GeometryTypes: Vec3
import PhysicalConstants.CODATA2018: c_0, m_e, e
using ..LaserTypes: TemporalProfiles, w, f

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

function Ex(z, r, par)
    @unpack E₀, w₀, k₀, z_R, φ₀ = par
    wz = w(z, par)

    E₀ * w₀/wz * exp(-im*k₀*z - r^2/wz^2 - im*((z*r^2)/(z_R*wz^2) + atan(z, z_R) + φ₀))
end

Ey(z, r, par) = par.ξy / par.ξx * Ex(z, r, par)

function Ez(x, y, z, r, par)
    @unpack k₀, z_R = par
    wz = w(z, par)

    2(im - z/z_R) / (k₀*wz^2) * (x*Ex(z, r, par) + y*Ey(z, r, par))
end

function E(x, y, z, t, par)
    r = hypot(x, y)

    E_x = Ex(z, r, par)
    E_y = Ey(z, r, par)
    E_z = Ez(x, y, z, r, par)

    real(Vec3(E_x, E_y, E_z) * g(z, t, par))
end

Bx(z, r, par) = -1/par.c * Ey(z, r, par)

By(z, r, par) = 1/par.c * Ex(z, r, par)

function Bz(x, y, z, r, par)
    @unpack k₀, z_R, c = par
    wz = w(z, par)

    2(im - z/z_R) / (c*k₀*wz^2) * (y*Ex(z, r, par) - x*Ey(z, r, par))
end

function B(x, y, z, t, par)
    r = hypot(x, y)

    B_x = Bx(z, r, par)
    B_y = By(z, r, par)
    B_z = Bz(x, y, z, r, par)

    real(Vec3(B_x, B_y, B_z) * g(z, t, par))
end

E(r, t, par) = E(r[1], r[2], r[3], t, par)
B(r, t, par) = B(r[1], r[2], r[3], t, par)

end  # module Gauss
