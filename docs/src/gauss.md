### [Gaussian profile](@id gaussian)

The `GaussLaser` type contains the parameters required for describing a
[Gaussian laser pulse](https://en.wikipedia.org/wiki/Gaussian_beam). The spatial profile of the electric field is given by
```math
\newcommand{\ii}{\mathrm{i}}                % imaginary unit
\begin{aligned}
    E_x (r,z) &= \xi_x E_0 \frac{w_0}{w} \exp\left[-\ii k z -\frac{r^2}{w^2} -\ii \left(\frac{k r^2}{2 R} - \arctg{\frac{z}{z_R}} - \phi_0\right)\right]\\
    E_y (r,z) &= \frac{\xi_y}{\xi_x} E_x(r,z) \\
    E_z (r,z) &= \frac{2 \left(\ii -\frac{z}{z_R}\right)}{kw^2 (z)} \left[xE_x (r,z) +yE_y(r,z)\right]\,,
\end{aligned}
```
where
- ``E_0`` is the amplitude of the electric field
- ``\xi_x`` and ``\xi_y`` give the polarization (choosing ``\xi_x=1`` and ``\xi_y=0`` for example, would give a linearly polarized field along the ``x`` axis while taking both ``1/\sqrt{2}`` would give circular polarization)
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
    B_z (r,z) &= \frac{2 \left(\ii -\frac{z}{z_R}\right)}{ckw^2 (z)} [yE_x (r,z) - xE_y(r,z)]
\end{aligned}
```
```@autodocs
Modules = [LaserTypes]
Pages = ["src/gauss.jl"]
Order = [:type]
```
!!! note Choice of normalization
The electric field is normalized such that in origin (with the default linear polarization) we obtain
```math
\vec{E}(0,0,0) = (E_0,0,0)\,,
```
when the temporal profile is ignored (or constant and equal to unity).