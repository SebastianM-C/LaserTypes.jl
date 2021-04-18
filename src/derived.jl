"""
    S(r, t, laser)

Compute the [Poynting vector](https://en.wikipedia.org/wiki/Poynting_vector) defined as
```math
    S(r, t) = \\frac{1}{μ₀} \\vec{E} × \\vec{B}
```
"""
function S(r, t, laser)
    ElectricField, MagneticField = EB(r, t, laser)

    1 / laser.μ₀ * (ElectricField × MagneticField)
end
