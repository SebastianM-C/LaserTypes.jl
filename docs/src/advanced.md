# Advanced usage

If you need more control over the parameters, you can use the constructors for the
laser types directly.

```@autodocs
Modules = [LaserTypes]
Pages = ["src/gauss.jl", "src/laguerre-gauss.jl"]
Order = [:type]
```
!!! note
    When specifying derived values, the corresponding independent parameters must also be given in order to avoid discrepancies. For example if you give `k`,
    you must also give the corresponding `λ`.