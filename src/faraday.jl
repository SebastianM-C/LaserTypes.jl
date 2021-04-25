# # 4-Potential

function EB(r, laser)
    x, y = r[1], r[2]
    coords = required_coords(laser, r)

    ElectricField = E(x, y, coords, laser)
    E_x = ElectricField[1]
    E_y = ElectricField[2]

    B_x = Bx(laser, E_y)
    B_y = By(laser, E_x)
    B_z = Bz(laser, coords, E_x, E_y, x, y)

    MagneticField = SVector{3}(B_x, B_y, B_z)

    return ElectricField, MagneticField
end

EB(r, t, laser) = real.(EB(r, laser) .* g(r[end], t, laser))

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
    c = laser.c
    r = x[begin+1:end]
    t = x[begin] / c
    (Ex, Ey, Ez), (Bx, By, Bz) = EB(r, t, laser)

    𝟘 = zero(Bx)
    return @SMatrix [𝟘     -Ex/c -Ey/c -Ez/c ;
                     Ex/c    𝟘    -Bz    By  ;
                     Ey/c    Bz    𝟘    -Bx  ;
                     Ez/c   -By    Bx    𝟘   ]
end
