# # Electric Field

function E(r, laser)
    coords = required_coords(laser, r)
    x, y = r[1], r[2]
    E(x, y, coords, laser)
end

E(r, t, laser) = real(E(r, laser) * g(r[3], t, laser))

function E(x, y, coords, laser)
    E_x = Ex(laser, coords)
    E_y = Ey(laser, E_x)
    E_z = Ez(laser, coords, E_x, E_y, x, y)

    SVector{3}(E_x, E_y, E_z)
end

Ey(laser, Ex) = laser.ξy / laser.ξx * Ex
