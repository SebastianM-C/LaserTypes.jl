# # Electric Field

function Ex(laser, x, y, z, r) end

Ey(laser, coords...) = laser.ξy / laser.ξx * Ex(laser, coords...)
Ey(laser, Ex) = laser.ξy / laser.ξx * Ex

function Ez(laser, x, y, z, r)
    E_x = Ex(laser, x, y, z, r)
    Ez(laser, E_x, Ey(laser, E_x), x, y, z, r)
end

function E(x, y, z, laser)
    r = hypot(x, y)

    E_x = Ex(laser, x, y, z, r)
    E_y = Ey(laser, E_x)
    E_z = Ez(laser, E_x, E_y, x, y, z, r)

    real(Vec3(E_x, E_y, E_z))
end

E(x, y, z, t, laser) = E(x, y, z, laser) * g(z, t, laser)

E(r, t, laser) = E(r[1], r[2], r[3], t, laser)
