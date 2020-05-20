# Laser initialization

Each supported laser has its own (Julia) type which contains all the parameters required to decribe the electromagnetic field configuration.
The simplest way to initialize a laser is via the `setup_laser` function.
```@autodocs
Modules = [LaserTypes]
Pages = ["src/setup.jl"]
```