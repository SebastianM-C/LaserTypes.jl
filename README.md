# LaserTypes
[build-img]: https://github.com/SebastianM-C/LaserTypes.jl/workflows/Run%20CI%20on%20master/badge.svg
[build-url]: https://github.com/SebastianM-C/LaserTypes.jl/actions
[codecov-img]: https://codecov.io/gh/SebastianM-C/LaserTypes.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/SebastianM-C/LaserTypes.jl
[docs-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-url]: http://SebastianM-C.github.io/LaserTypes.jl/dev/

[![build status][build-img]][build-url]
[![codecov][codecov-img]][codecov-url]
[![docs][docs-img]][docs-url]

This package aims to provide a common interface for different laser types. For the moment only the Gaussian pulse is supported, but in the future, Laguerre-Gauss and Bessel beams will be added.

Each laser type has its own type. For example for the Gaussian laser pusle use:
```julia
using LaserTypes
using Unitful
using StaticArrays

p = GaussLaser()
```
This will give the functions for the values of the electromagnetic field at a space-time point specifed by `r,t`. For example, to evaluate
the electric field at the origin use
```julia
x₀ = SVector{3}(0u"μm",0u"μm",0u"μm")
t₀ = 0u"s"

E(x₀, t₀, p)
```

The package does not impose the use of units, but by default the laser parameters are initialized in SI.
For example with regular `Float64` number, the laser initialization would look like this:
```julia
c = 137.035
q = -1
m = 1
λ = 15117.8089
w0 = 944863.062
τ = 744.144
t₀ = 0
p = GaussLaser(c=c, q=q, m_q=m, λ=λ, w₀=w0, profile=GaussianProfile(c=c,τ=τ))
```
We can vizualize the intensity of the created electric field with Makie.jl like this:
```julia
using Makie
using LinearAlgebra

f(x,y) = norm(E(Point3f0(x*10^6,y*10^6,p.z_F), 1, p))
surface(-5:0.1:5, -5:0.1:5, f)
```
![gauss](assets/gauss.png)
