# # Electric Field

function E(r, laser)
    coords = required_coords(laser, r)
    x, y = r[1], r[2]
    E(x, y, coords, laser)
end

E(r, t, laser) = E(r, laser) * real(g(r[3], t, laser))

function E(x, y, coords, laser)
    E_x = Ex(laser, coords)
    E_y = Ey(laser, E_x)
    E_z = Ez(laser, coords, E_x, E_y, x, y)

    real(Vec3(E_x, E_y, E_z))
end

E(x, y, z, t, laser) = E(Vec3(x, y, z), t, laser)

Ey(laser, Ex) = laser.ξy / laser.ξx * Ex
