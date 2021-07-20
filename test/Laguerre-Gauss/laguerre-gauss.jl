using SafeTestsets, Test

@safetestset "Values at origin" begin include("origin.jl") end
@safetestset "Orientation" begin include("orientation.jl") end
@safetestset "Field values compared with symbolic results" begin include("field_values.jl") end
@safetestset "Divergence in Polynomial Roots" begin include("roots.jl") end
@safetestset "LG(p = 0, m = 0, x, y, z, t) == G(x, y, z, t)" begin include("gauss.jl") end
@safetestset "Comparison to Symbolic Values" begin include("field_values.jl") end
