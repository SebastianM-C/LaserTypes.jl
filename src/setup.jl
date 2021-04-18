"""
    setup_laser(laser, units; τ=nothing, kwargs...)

Initialize the specified `laser` using the default parameteres in the specified
`units`. Supported units are
- `:SI`: Values are in SI, the numeric types are regular numbers (`Float64`)
- `:SI_unitful`: Values are in SI, but using [Unitful.jl](https://github.com/PainterQubits/Unitful.jl)
- `:atomic`: Values are in atomic units, the numeric types are regular numbers (`Float64`)
- `:atomic_unitful`: Values are in atomic units, but using [UnitfulAtomic.jl](https://github.com/sostock/UnitfulAtomic.jl)

The keyword arguments can be used to give specific values to the parameters instead
of using the defaults.
You can specify parameteres such as the wavelength and beam waist via `λ` and `w₀`.
The duration of the pulse (assuming Gaussian temporal profile) can be given via `τ`.
Be sure to use the same units as the ones provided via `units`.
You can also specify dimensionless parameteres such as the normalized vector potential (`a₀`),
the initial phase (`ϕ₀`) and the polarization (`ξx` and `ξy`).
See the docomuntation for each laser type for more details on the supported arguments.
"""
function setup_laser(laser, units; τ=nothing, kwargs...)
    _c, _q, _m_q, _μ₀ = fundamental_constants(Val(units))
    _λ, _w₀ = common_parameters(Val(units))
    c = get(kwargs, :c, _c)
    q = get(kwargs, :q, _q)
    m_q = get(kwargs, :m_q, _m_q)
    μ₀ = get(kwargs, :m_q, _μ₀)
    λ = get(kwargs, :λ, _λ)
    w₀ = get(kwargs, :w₀, _w₀)

    if τ ≡ nothing
        τ = duration(Val(units))
    end
    profile = get(kwargs, :profile, GaussProfile(c=c, τ=τ))

    others = Dict{Symbol,Any}()
    excluded = [:c, :q, :m_q, :μ₀, :λ, :w₀, :profile]
    for (k,v) in kwargs
        if k ∉ excluded
            push!(others, k=>v)
        end
    end

    laser(;c=c, q=q, m_q=m_q, μ₀=μ₀, λ=λ, w₀=w₀, profile=profile, others...)
end

fundamental_constants(::Val{:SI_unitful}) = c_0, -e, m_e, μ_0
fundamental_constants(::Val{:SI}) = ustrip.(fundamental_constants(Val(:SI_unitful)))
fundamental_constants(::Val{:atomic_unitful}) = auconvert.(fundamental_constants(Val(:SI_unitful)))
fundamental_constants(::Val{:atomic}) = ustrip.(fundamental_constants(Val(:atomic_unitful)))

common_parameters(::Val{:SI_unitful}) = 0.8u"μm", 58.0u"μm"
common_parameters(::Val{:SI}) = ustrip.(u"m", common_parameters(Val(:SI_unitful)))
common_parameters(::Val{:atomic_unitful}) = auconvert.(common_parameters(Val(:SI_unitful)))
common_parameters(::Val{:atomic}) = ustrip.(common_parameters(Val(:atomic_unitful)))

duration(::Val{:SI_unitful}) = 18.0u"fs"
duration(::Val{:SI}) = ustrip.(u"s", duration(Val(:SI_unitful)))
duration(::Val{:atomic_unitful}) = auconvert.(duration(Val(:SI_unitful)))
duration(::Val{:atomic}) = ustrip.(duration(Val(:atomic_unitful)))
