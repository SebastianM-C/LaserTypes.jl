module LaserTypes

export E, B, Gauss, LaguerreGauss, TemporalProfiles

using Parameters

function w(z, par)
    @unpack w₀, z_R = par

    w₀ * √(1 + (z/z_R)^2)
end

R(z, z_R) = z + z_R^2 / z

"""
    g(z, t, par)

The time dependent part of the laser pulse.
"""
function g(z, t, par)
    @unpack envelope, ω = par

    exp(im*ω*t) * envelope(z, t, par)
end

E(x, y, z, t, par) = E(x, y, z, par) * real(g(z, t, par))
B(x, y, z, t, par) = B(x, y, z, par) * real(g(z, t, par))

E(r, t, par) = E(r[1], r[2], r[3], t, par)
B(r, t, par) = B(r[1], r[2], r[3], t, par)

function EB(r, t, par)
    x, y, z = r[1], r[2], r[3]
    ElectricField = E(x, y, z, par)
    Ex = ElectricField[1]
    Ey = ElectricField[2]

    B_x = Bx(Ey, par)
    B_y = By(Ex, par)
    B_z = Bz(Ex, Ey, x, y, z, par)

    MagneticField = real(Vec3(B_x, B_y, B_z))

    temporal_part = real(g(z, t, par))
    ElectricField *= temporal_part
    MagneticField *= temporal_part

    return ElectricField, MagneticField
end

include("envelopes.jl")
include("gauss.jl")
include("laguerre-gauss.jl")

end # module
