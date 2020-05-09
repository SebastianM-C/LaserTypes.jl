### Laguerre-Gauss profile

The `LaguerreGaussLaser` type contains the parameters required for describing a
Laguerre-Gauss laser pulse. The spatial profile of the electric field is given by
```math
\newcommand{\ii}{\mathrm{i}}                % imaginary unit
\begin{aligned}
    E
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
