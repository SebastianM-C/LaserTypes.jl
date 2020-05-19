using LaserTypes
using Test, SafeTestsets

@testset "LaserTypes.jl" begin
    @safetestset "Setup" begin include("setup.jl") end
    @safetestset "Dimensions" begin include("dimensions.jl") end
    @safetestset "Gauss" begin include("gauss.jl") end
    @safetestset "Laguerre-Gauss" begin include("laguerre-gauss.jl") end
end
