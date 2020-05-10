# Temporal Profiles

The complete expression of the electromaginetic fields describing the laser pulse
is given by the combination of the spatial part together with the time dependence.
```@autodocs
Modules = [LaserTypes]
Pages = ["LaserTypes.jl"]
```

The `TemporalProfiles` profiles module contains some predefined envelopes.

### Constant profile

This is the trivial profile
```math
envelope(z, t) = 1
```
which gives an infinite duration for the electromagnetic field configuration (we cannot
call this a laser pusle since a pulse implicitely has a finite duration).

### Gaussian profile

This envelope provides a finite duration for the laser pulse and thus can provide a more
realistic description of an actual laser pulse.
```math
envelope(z, t) = \exp\left[-\left(\frac{t-\frac{z-z_F}{c}}{\tau_0}\right)^2\right]\,,
```
where
- `τ₀` is the duration of the pulse (FWHM) and has the default value 18.02fs
- `z_F` is the initial position of the intensity peak and has the default value `-4*τ₀*c`

### Quasi-rectangular profile

This envelope
