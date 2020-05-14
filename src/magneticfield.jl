# # Magnetic Field

function B(r, laser)
    coords = required_coords(laser, r)
    x, y = r[1], r[2]
    B(x, y, coords, laser)
end

B(r, t, laser) = B(r, laser) * g(r[3], t, laser)

function B(x, y, coords, laser)
    E_x = Ex(laser, coords)
    E_y = Ey(laser, E_x)

    B_x = Bx(laser, E_y)
    B_y = By(laser, E_x)
    B_z = Bz(laser, coords, E_x, E_y, x, y)

    real(Vec3(B_x, B_y, B_z))
end

B(x, y, z, t, laser) = B(Vec3(x, y, z), t, laser)

Bx(laser, Ey) = -1/laser.c * Ey
By(laser, Ex) = 1/laser.c * Ex
