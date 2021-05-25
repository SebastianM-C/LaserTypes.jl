# # Magnetic Field

function B(r, t, laser)
    R = geometry(laser).rotation_matrix
    r′ = rotate_coords(R, r)
    inv_c = immutable_cache(laser, :inv_c)
    if unit(eltype(r)) ≠ NoUnits
        λ = laser.λ
        r′ = uconvert.(unit(λ), r′)
    end

    MagneticField = real(B(r′, laser) * g(r′[3], t, laser; inv_c))

    return inv_rotate(R, MagneticField)
end

function B(x::SVector{4}, laser::AbstractLaser)
    inv_c = immutable_cache(laser, :inv_c)
    r = x[SVector{3}(2,3,4)]
    t = x[begin] * inv_c
    B(r, t, laser)
end

function B(r, laser::AbstractLaser)
    @assert length(r) == 3 "The laser is only defined in 3D"
    fill!(laser.cache, r)
    coords = required_coords(laser, r)

    E_x = Ex(laser, coords)
    update_cache!(laser, :Ex, E_x)
    E_y = Ey(laser, coords)
    update_cache!(laser, :Ey, E_y)

    B_x = Bx(laser, coords)
    B_y = By(laser, coords)
    B_z = Bz(laser, coords)

    SVector{3}(B_x, B_y, B_z)
end

@inline function Bx(laser, coords)
    c = fundamental_constants(laser, :c)
    Ey = mutable_cache(laser).Ey

    -1/c * Ey
end

@inline function By(laser, coords)
    c = fundamental_constants(laser, :c)
    Ex = mutable_cache(laser).Ex

    1/c * Ex
end
