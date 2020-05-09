# # Gaussian profile

@with_kw struct GaussLaser{V,Q,M,L,F,C,T,P,W,K,E}
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

@inline Ex(laser::GaussLaser, x, y, z, r) = Ex(laser, z, r)

function Ex(laser::GaussLaser, z, r)
    @unpack E₀, w₀, k, z_R, φ₀, ξx = laser
    wz = w(z, laser)
    Rz = R(z, z_R)

    ξx * E₀ * w₀/wz * exp(-im*k*z - (r/wz)^2 - im*((k*r^2)/(2Rz) - atan(z, z_R) - φ₀))
end

@inline Ez(laser::GaussLaser, Ex, Ey, x, y, z, r) = Ez(laser, Ex, Ey, x, y, z)

function Ez(laser::GaussLaser, Ex, Ey, x, y, z)
    @unpack k, z_R = laser
    wz = w(z, laser)

    2(im - z/z_R) / (k*wz^2) * (x*Ex + y*Ey)
end

@inline Bz(laser::GaussLaser, Ex, Ey, x, y, z, r) = Bz(laser, Ex, Ey, x, y, z)

function Bz(laser::GaussLaser, Ex, Ey, x, y, z)
    @unpack k, z_R, c = laser
    wz = w(z, laser)

    2(im - z/z_R) / (c*k*wz^2) * (y*Ex - x*Ey)
end
