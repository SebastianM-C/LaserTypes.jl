# Temporal Profiles

The complete expression of the electromaginetic fields describing the laser pulse
is given by the combination of the spatial part together with the time dependence.
```@docs
LaserTypes.g(::Any, ::Any, ::Any)
```

There are several types that provide profiles to the temporal envelope:
```@autodocs
Modules = [LaserTypes]
Pages = ["src/envelopes.jl"]
Order = [:type]
```
