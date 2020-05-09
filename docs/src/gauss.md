### Gaussian profile

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

The `GaussLaser` is defined by the following independent parameters
- `c` is the speed of light in vaccum, with the default value being in SI (`c_0` from the CODATA2018 in the [PhysicalConstants](https://github.com/JuliaPhysics/PhysicalConstants.jl) package)
- `q` is the electric charge, with the default value being the one for the electron in SI (`-e` from the CODATA2018 in the [PhysicalConstants](https://github.com/JuliaPhysics/PhysicalConstants.jl) package)
- `m_q` is the mass of the charge, with the default value being the one for the electron in SI (`m_e` from the CODATA2018 in the [PhysicalConstants](https://github.com/JuliaPhysics/PhysicalConstants.jl) package)
- `λ` is the laser wavelangth with the default value 0.8μm
- `a₀` is the normalized vector potential (defined as ``a_0=\frac{eA}{m_e c^2}``)
- `φ₀` is the initial phase with the default value 0.0
- `w₀` is the beam radius at the Rayleigh range or [beam waist](https://en.wikipedia.org/wiki/Gaussian_beam#Beam_waist) with the default value 58.0μm
- `ξx` and `ξy` give the polarization and have the default value `1.0 + 0im` and `0.0 + 0im`
- `τ₀` is the ... and has the default value 18.02fs
- `z_F` is the ... and has the default value `-4*τ₀*c`
- `envelope` is the temporal envelope of the pulse and the default one is a Gaussian one

During the initialization of the `GaussLaser` type, some useful derived values are
also computed
- `ω` is the angular frequency and is given by ``2\pi \frac{c}{\lambda}``
- `k` is the wavenumber and is given by ``\frac{2\pi}{\lambda}``
- `z_R` is the [Rayleigh range](https://en.wikipedia.org/wiki/Rayleigh_length) and is given by ``\frac{k w_0^2}{2}``
- `T₀` is the laser period and is given by ``\frac{2\pi}{\omega}``
- `E₀` is the amplitude of the electric field and is given by ``a_0\frac{m_q c \omega}{q}``

The electric field is normalized such that in origin (with the default linear polarization) we obtain
```math
\vec{E}(0,0,0) = (E_0,0,0)\,,
```
when the temporal profile is ignored (or constant and equal to unity).
