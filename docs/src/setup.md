# Laser initialization

Each supported laser has its own (Julia) type which contains all the parameters required to decribe the electromagnetic field configuration.
The simplest way to initialize a laser is via the `setup_laser` function.
```@autodocs
Modules = [LaserTypes]
Pages = ["src/setup.jl"]
```

For example, you can create a Gaussian laser with ``a_0 = 2``,
``λ = 800 \text{nm}`` and waist size ``w_0 = 50 \text{μm}``
and a duration ``τ = 18 \text{fs}`` in the following ways

```@example setup
using LaserTypes
using Unitful

λ = 800u"nm"
w0 = 50u"μm"
τ = 18.0u"fs"

laser = setup_laser(GaussLaser, :SI_unitful; λ, w₀=w0, a₀=2.0, τ)
```

or

```@example setup
import PhysicalConstants.CODATA2018: c_0

ω = 2π * c_0 / λ

laser = setup_laser(GaussLaser, :SI_unitful; ω, w₀=w0, a₀=2.0, τ)
```

or

```@example setup
k = 2π / λ

laser = setup_laser(GaussLaser, :SI_unitful; k, w₀=w0, a₀=2.0, profile = GaussProfile, τ)
```

or

```@example setup
z_R = w0^2 * k / 2

laser = setup_laser(GaussLaser, :SI_unitful; k, z_R, a₀=2.0, profile = GaussProfile(;τ, z₀=0u"nm"))
```
