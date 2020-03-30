# LaserTypes
[![build status](https://github.com/SebastianM-C/LaserTypes.jl/workflows/Run%20CI%20on%20master/badge.svg)](https://github.com/SebastianM-C/LaserTypes.jl/actions)
[![codecov](https://codecov.io/gh/SebastianM-C/LaserTypes.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SebastianM-C/LaserTypes.jl)

This package aims to provide a common interface for different laser types. For the moment only the Gaussian pulse is supported, but in the future, Laguerre-Gauss and Bessel beams will be added.

Each laser type has its own module. For example for the Gaussian laser pusle use:
```julia
using LaserTypes
using Unitful
using StaticArrays

p = Gauss.LaserParams()

E = Gauss.E
B = Gauss.B
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
τ0 = 744.144
t₀ = 0
p = Gauss.LaserParams(c=c, q=q, m_q=m, λ=λ, w₀=w0, τ₀=τ0)
```
We can vizualize the intensity of the created electric field with Makie.jl like this:
```julia
f(x,y) = norm(E(Point3f0(x*10^6,y*10^6,p.z_F), 1, p))
surface(-5:0.1:5, -5:0.1:5, f)
```
