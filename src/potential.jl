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
