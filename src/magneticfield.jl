# # Magnetic Field

Bx(laser, c1, c2, rest_of_coords...) = -1/laser.c * Ey(laser, c1, c2, rest_of_coords...)
Bx(laser, Ey) = -1/laser.c * Ey

By(laser, c1, c2, rest_of_coords...) = 1/laser.c * Ex(laser, c1, c2, rest_of_coords...)
By(laser, Ex) = 1/laser.c * Ex

function Bz(laser, x, y, z, r)
    E_x = Ex(laser, x, y, z, r)
    Bz(laser, E_x, Ey(laser, E_x), x, y, z, r)
end

function B(x, y, z, laser)
    r = hypot(x, y)

    E_x = Ex(laser, x, y, z, r)
    E_y = Ey(laser, E_x)

    B_x = Bx(laser, E_y)
    B_y = By(laser, E_x)
    B_z = Bz(laser, E_x, E_y, x, y, z, r)

    real(Vec3(B_x, B_y, B_z))
end

B(x, y, z, t, laser) = B(x, y, z, laser) * real(g(z, t, laser))

B(r, t, laser) = B(r[1], r[2], r[3], t, laser)
