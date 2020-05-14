### Laguerre-Gauss profile

The `LaguerreGaussLaser` type contains the parameters required for describing a
[Laguerre-Gauss laser](https://en.wikipedia.org/wiki/Gaussian_beam#Laguerre-Gaussian_modes) pulse. The spatial profile of the electric field is given by
```math
\newcommand{\ii}{\mathrm{i}}                % imaginary unit
\begin{aligned}
    E_x(r,z,φ) &= ξ_x E_g N_{pm} \left(\frac{\sqrt{2}r}{w(z)}\right)^{|m|} ₁F₁\left(-p, |m|+1; 2\left(\frac{r}{w(z)}\right)^2\right) \exp\left[\ii\left((2p+|m|)\arctg{\frac{z}{z_R}}-mφ - \phi_0\right)\right] \\
    E_y(r,z,φ) &= \frac{ξ_y}{ξ_x} E_x(r,z,φ) \\
    E_z(r,z,φ) &= -\frac{\ii}{k} \left\{\left[-2\frac{1+\ii\frac{z}{z_0}}{w^2(z)} + \frac{4p}{(|m|+1)w(z)^2}\ ₁F₁^{-1}\left(-p, |m|+1; 2\left(\frac{r}{w(z)}\right)^2\right)\right](x E_x + y E_y) - \frac{|m|}{x+\ii y} (E_x \mp \ii E_y)\right\}
\end{aligned}
```
where
- ``E_g`` is the amplitude of the electric field form the [Gaussian profile](@ref gaussian), ``E_x(r,z)``
- ``N_{pm}`` is a normalization factor given by ``\sqrt{(p+1)_{|m|}}``, with ``(x)_n`` the [Pochhammer symbol](https://mathworld.wolfram.com/PochhammerSymbol.html)
- ``₁F₁(α,β;z)`` is the [confluent hypergeometric function of the first kind](https://mathworld.wolfram.com/ConfluentHypergeometricFunctionoftheFirstKind.html)
- ``p`` is the radial index of the mode, ``p ∈ ℤ, p ≥ 0``
- ``m`` is the azimuthal index of the mode, ``m ∈ ℤ``
- ``ξ_x`` and ``ξ_y`` give the polarization (choosing ``ξ_x=1`` and ``ξ_y=0`` for example, would give a linearly polarized field along the ``x`` axis while taking them ``1/\sqrt{2}`` and ``±\mathrm{i}/\sqrt{2}`` would give right and left-handed circular polarization)
- ``w = w(z)`` is the beam radius at ``z``
- ``w_0`` is the [beam waist](https://en.wikipedia.org/wiki/Gaussian_beam#Beam_waist)
- ``k`` is the wavenumber
- ``R`` is the radius of curvature
- ``\arctg{\frac{z}{z_R}}`` is the Gouy phase
- ``\phi_0`` is the initial phase

The magnetic field field is given by
```math
\newcommand{\ii}{\mathrm{i}}                % imaginary unit
\begin{aligned}
    B_x(r,z) &= -\frac{1}{c}E_y(r,z)\\
    B_y(r,z) &= \frac{1}{c}E_x(r,z)\\
    B_z (r,z) &= -\frac{\ii}{ck} \left\{\left[-2\frac{1+\ii\frac{z}{z_0}}{w^2(z)} + \frac{4p}{(|m|+1)w(z)^2}\ ₁F₁^{-1}\left(-p, |m|+1; 2\left(\frac{r}{w(z)}\right)^2\right)\right](y E_x + x E_y) - \frac{|m|}{x+\ii y} (E_y \mp \ii E_x)\right\}
\end{aligned}
```
```@autodocs
Modules = [LaserTypes]
Pages = ["src/laguerre-gauss.jl"]
```
!!! note "Choice of normalization"

    The normalization was chosen such that the normalization factor ``N_{pm}``
    is adimensional and the result is consistent with the choice for the Gaussian
    mode. This way it is easy to obtain the electric field from the mode, as the
    dimensions are given through ``E_0``.

    The Laguerre-Gauss modes are expressed via the Gauss ones in order to reuse
    code and together with the choice for ``N_{pm}`` ensures consistent normalization
    preferences.

!!! note

    In the expression of the Laguerre-Gauss modes above, the [confluent hypergeometric function of the first kind](https://mathworld.wolfram.com/ConfluentHypergeometricFunctionoftheFirstKind.html) was used instead of the
    standard [associated Laguerre polynomials](https://en.wikipedia.org/wiki/Laguerre_polynomials#Generalized_Laguerre_polynomials). This helps reduce the number of
    factorials required for the normalization constant. The [relation](https://en.wikipedia.org/wiki/Laguerre_polynomials#Relation_to_hypergeometric_functions)
    between the two is given by
    ```math
    L_n^α(z) = \frac{(α+1)_n}{n!} ₁F₁(-n,α+1;z),
    ```
    where ``(x)_n`` is the [Pochhammer symbol](https://mathworld.wolfram.com/PochhammerSymbol.html)
