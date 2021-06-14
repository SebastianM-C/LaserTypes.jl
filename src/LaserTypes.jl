# # LaserTypes

module LaserTypes

export E, B,
    GaussLaser, LaguerreGaussLaser,
    ConstantProfile, GaussProfile, Cos²Profile, QuasiRectangularProfile,
    setup_laser, Fμν, S

using Unitful
using UnitfulAtomic
using UnPack
using LinearAlgebra
using HypergeometricFunctions
using StaticArrays
using CoordinateTransformations
using AutoHashEquals
using ParallelProcessingTools
import PhysicalConstants.CODATA2018: c_0, e, m_e, ε_0, μ_0

const _₁F₁ = HypergeometricFunctions.drummond1F1
const pochhammer = HypergeometricFunctions.pochhammer

abstract type AbstractLaser end
abstract type AbstractTemporalProfile end

include("envelopes.jl")
include("constants.jl")
include("coords.jl")
include("polarization.jl")
include("electricfield.jl")
include("magneticfield.jl")
include("faraday.jl")
include("gauss.jl")
include("laguerre-gauss.jl")
include("setup.jl")
include("derived.jl")
include("utils.jl")
include("show.jl")

end # module
