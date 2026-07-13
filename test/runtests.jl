using LaserTypes
using Test, SafeTestsets

@testset "LaserTypes.jl" begin
    @safetestset "Setup" begin include("setup.jl") end
    @safetestset "Dimensions" begin include("dimensions.jl") end
    @safetestset "Gauss" begin include("gauss.jl") end
    @safetestset "Laguerre-Gauss" begin include("Laguerre-Gauss/laguerre-gauss.jl") end
    @safetestset "Fμν" begin include("faraday.jl") end
    @safetestset "Derived" begin include("derived.jl") end
    @safetestset "Threads" begin include("threads.jl") end
    @safetestset "Show" begin include("show.jl") end
    # JET reports upstream runtime dispatch on Julia pre-releases (e.g. Base
    # printing internals and StaticArrays construction on 1.13-DEV), so only
    # run the QA tests on release versions.
    if isempty(VERSION.prerelease)
        @safetestset "QA" begin include("qa.jl") end
    else
        @info "Skipping JET QA tests on pre-release Julia" VERSION
    end
end
