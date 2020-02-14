using LaserTypes
using Test, SafeTestsets

@testset "LaserTypes.jl" begin
    @safetestset "Gauss" begin include("gauss.jl") end
end
