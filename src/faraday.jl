# # 4-Potential

function EB(r, t, laser)
    R = geometry(laser).rotation_matrix
    r′ = rotate_coords(R, r)
    inv_c = immutable_cache(laser, :inv_c)
    if unit(eltype(r)) ≠ NoUnits
        λ = laser.λ
        r′ = uconvert.(unit(λ), r′)
    end

    E_B = real.(EB(r′, laser) .* g(r′[end], t, laser; inv_c))

    return inv_rotate.((R,), E_B)
end

function EB(r, laser)
    @assert length(r) == 3 "The laser is only defined in 3D"
    fill!(mutable_cache(laser), r)
    coords = required_coords(laser, r)

    E_x = Ex(laser, coords)
    update_cache!(laser, :Ex, E_x)
    E_y = Ey(laser, coords)
    update_cache!(laser, :Ey, E_y)
    E_z = Ez(laser, coords)
    update_cache!(laser, :Ez, E_z)

    ElectricField = SVector{3}(E_x, E_y, E_z)

    B_x = Bx(laser, coords)
    B_y = By(laser, coords)
    B_z = Bz(laser, coords)

    MagneticField = SVector{3}(B_x, B_y, B_z)

    return ElectricField, MagneticField
end

"""
    Fμν(x, laser)

Compute [the electromagnetic tensor](https://en.wikipedia.org/wiki/Electromagnetic_tensor) ``F`` in
covarinat matrix form at the spacetime 4-vector specified by `x`.
The ``x = (x^0 ≡ ct, x^1, x^2, x^3)`` convention is used.
```math
F^{μν} =
\\left[\\begin{array}{cccc}
    0   &  -E_x / c  & -E_y / c & -E_z / c \\\\
E_x / c &       0    &   -B_z   &    B_y   \\\\
E_y / c &      B_z   &     0    &   -B_x   \\\\
E_z / c &     -B_y   &    B_z   &     0
\\end{array}\\right]
```
"""
function Fμν(x, laser)
    c⁻¹ = immutable_cache(laser, :inv_c)
    r = x[begin+1:end]
    t = x[begin] * c⁻¹
    (Ex, Ey, Ez), (Bx, By, Bz) = EB(r, t, laser)

    𝟘 = zero(Bx)
    return @SMatrix [𝟘      -Ex*c⁻¹  -Ey*c⁻¹ -Ez*c⁻¹ ;
                     Ex*c⁻¹   𝟘      -Bz      By     ;
                     Ey*c⁻¹   Bz      𝟘      -Bx     ;
                     Ez*c⁻¹  -By      Bx      𝟘      ]
end
