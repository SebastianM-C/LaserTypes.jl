module Gauss

using Unitful
using Parameters
using GeometryTypes: Vec3
import PhysicalConstants.CODATA2018: c_0, m_e, e
using ..LaserTypes: LaserTypes, TemporalProfiles, w, g, R

@with_kw struct GaussParams{V,Q,M,L,F,C,T,P,W,K,E}
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
    # derived values
    ω::W = 2π * c / λ; @assert ω ≈ 2π * c / λ
    k::K = 2π / λ; @assert k ≈ 2π / λ
    z_R::L = w₀^2 * k / 2; @assert z_R ≈ w₀^2 * k / 2
    T₀::T = uconvert(unit(τ₀), 2π / ω); @assert T₀ ≈ 2π / ω
    E₀::E = a₀ * m_q * c * ω / abs(q); @assert E₀ ≈ a₀ * m_q * c * ω / abs(q)
end

function Ex(z, r, par)
    @unpack E₀, w₀, k, z_R, φ₀, ξx = par
    wz = w(z, par)
    Rz = R(z, z_R)

    ξx * E₀ * w₀/wz * exp(-im*k*z - (r/wz)^2 - im*((k*r^2)/(2Rz) - atan(z, z_R) - φ₀))
end

Ey(z, r, par) = par.ξy / par.ξx * Ex(z, r, par)
Ey(Ex, par) = par.ξy / par.ξx * Ex

function Ez(x, y, z, r, par)
    E_x = Ex(z, r, par)
    Ez(E_x, Ey(E_x, par), x, y, z, par)
end

function Ez(Ex, Ey, x, y, z, par)
    @unpack k, z_R = par
    wz = w(z, par)

    2(im - z/z_R) / (k*wz^2) * (x*Ex + y*Ey)
end

function LaserTypes.E(x, y, z, par::GaussParams)
    r = hypot(x, y)

    E_x = Ex(z, r, par)
    E_y = Ey(E_x, par)
    E_z = Ez(E_x, E_y, x, y, z, par)

    real(Vec3(E_x, E_y, E_z))
end

Bx(z, r, par) = -1/par.c * Ey(z, r, par)
Bx(Ey, par) = -1/par.c * Ey

By(z, r, par) = 1/par.c * Ex(z, r, par)
By(Ex, par) = 1/par.c * Ex

function Bz(x, y, z, r, par)
    @unpack k, z_R, c = par
    wz = w(z, par)

    2(im - z/z_R) / (c*k*wz^2) * (y*Ex(z, r, par) - x*Ey(z, r, par))
end

function Bz(Ex, Ey, x, y, z, par)
    @unpack k, z_R, c = par
    wz = w(z, par)

    2(im - z/z_R) / (c*k*wz^2) * (y*Ex - x*Ey)
end

function LaserTypes.B(x, y, z, par::GaussParams)
    r = hypot(x, y)

    E_x = Ex(z, r, par)
    E_y = Ey(E_x, par)

    B_x = Bx(E_y, par)
    B_y = By(E_x, par)
    B_z = Bz(E_x, E_y, x, y, z, par)

    real(Vec3(B_x, B_y, B_z))
end

end  # module Gauss
