# # Electric Field

function E(r, t, laser)
    R = geometry(laser).rotation_matrix
    r′ = rotate_coords(R, r)
    inv_c = immutable_cache(laser, :inv_c)
    if unit(eltype(r)) ≠ NoUnits
        λ = laser.λ
        r′ = uconvert.(unit(λ), r′)
    end

    ElectricField = real(E(r′, laser) * g(r′[3], t, laser; inv_c))

    return inv_rotate(R, ElectricField)
end

function E(x::SVector{4}, laser::AbstractLaser)
    inv_c = immutable_cache(laser, :inv_c)
    r = x[SVector{3}(2,3,4)]
    t = x[begin] * inv_c
    E(r, t, laser)
end

function E(r, laser::AbstractLaser)
    @assert length(r) == 3 "The laser is only defined in 3D"
    fill!(mutable_cache(laser), r)
    coords = required_coords(laser, r)

    E_x = Ex(laser, coords)
    update_cache!(laser, :Ex, E_x)
    E_y = Ey(laser, coords)
    update_cache!(laser, :Ey, E_y)
    E_z = Ez(laser, coords)
    update_cache!(laser, :Ez, E_z)

    return SVector{3}(E_x, E_y, E_z)
end

@inline function Ey(laser, coords)
    @unpack ξx, ξy = polarization(laser)
    @unpack Ex = mutable_cache(laser)

    return ξy / ξx * Ex
end
