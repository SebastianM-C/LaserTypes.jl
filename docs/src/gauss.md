### Gaussian profile

The `GaussLaser` type contains the parameters required for describing a
Gaussian laser pulse. The spatial profile of the electric field is given by
```math
\newcommand{\ii}{\mathrm{i}}                % imaginary unit
\begin{aligned}
    E_x (z,r) &= \xi_x E_0 \frac{w_0}{w} \exp\left(-\frac{r^2}{w^2}  -\ii k z -\ii \frac{k r^2}{2 R} + \ii \phi_0\right)\\
    E_y (z,r) &= \xi_y E_0 \frac{w_0}{w} \exp\left(-\frac{r^2}{w^2}  -\ii k z -\ii \frac{k r^2}{2 R} + \ii \phi_0\right)\\
    E_z (z,r) &= \frac{2 \left(\ii -\frac{z}{z_0}\right)}{kw^2 (z)} \left[xE_x (r,z) +yE_y(r,z)\right]\,,
\end{aligned}
```
where
- ``E_0`` is the amplitude of the electric field
- ``\xi_x`` and ``\xi_y`` give the polarization (choosing ``\xi_x=1`` and ``\xi_y=0`` for example, would give a linearly polarized field along the ``x`` axis while taking both ``1/\sqrt{2}`` would give circular polarization)
- ``w = w(z)`` is the beam radius at ``z``
- ``w_0`` is the beam radius at the Rayleigh range ``z=z_R``
- ``k`` is the wavenumber
- ``R`` is the radius of curvature
- ``\phi_0`` is the Gouy phase

The magnetic field
