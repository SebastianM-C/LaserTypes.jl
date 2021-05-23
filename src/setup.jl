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
function setup_laser(laser, units; kwargs...)
    _λ, _w₀ = common_parameters(Val(units))
    constants = FundamentalConstants(units)
    c = constants.c

    if haskey(kwargs, :λ)
        λ = kwargs[:λ]
    elseif haskey(kwargs, :ω)
        ω = kwargs[:ω]
        λ = 2π*c / ω
    elseif haskey(kwargs, :k)
        k = kwargs[:k]
        λ = 2π/k
    else
        λ = _λ
    end

    if haskey(kwargs, :w₀)
        w₀ = kwargs[:w₀]
    elseif haskey(kwargs, :z_R)
        z_R = kwargs[:z_R]
        k = 2π / λ
        w₀ = √(2z_R / k)
    else
        w₀ = _w₀
    end

    if haskey(kwargs, :a₀)
        a₀ = kwargs[:a₀]
    elseif haskey(kwargs, :E₀)
        E₀ = kwargs[:E₀]
        q = abs(constants.q)
        mₑ = constants.mₑ
        ω = 2π * c / λ
        a₀ = E₀ * q / (mₑ * c * ω)
    else
        a₀ = 1
    end

    τ = get(kwargs, :τ, duration(Val(units)))
    t₀ = get(kwargs, :t₀, zero(τ))
    z₀ = get(kwargs, :z₀, zero(λ))

    @debug "profile"
    if haskey(kwargs, :profile)
        prof = kwargs[:profile]
        if prof isa Type
            if prof isa Type{<:GaussProfile} || prof isa Type{<:Cos²Profile}
                profile = prof(;τ, t₀, z₀)
            elseif prof isa Type{<:QuasiRectangularProfile}
                Δt = get(kwargs, :Δt, 10τ)
                profile = prof(;τ, t₀, z₀, Δt)
            else
                profile = prof()
            end
        else
            profile = prof
        end
    else
        @debug "default profile"
        profile = GaussProfile(;τ, z₀)
    end

    others = Dict{Symbol,Any}()
    processed = [:λ, :ω, :k, :w₀, :z_R, :a₀, :E₀, :τ, :t₀, :z₀, :profile]
    for (k,v) in kwargs
        if k ∉ processed
            push!(others, k=>v)
        end
    end

    laser(units; λ, w₀, a₀, profile, others...)
end

common_parameters(::Val{:SI_unitful}) = 0.8u"μm", 58.0u"μm"
common_parameters(::Val{:SI}) = ustrip.(u"m", common_parameters(Val(:SI_unitful)))
common_parameters(::Val{:atomic_unitful}) = auconvert.(common_parameters(Val(:SI_unitful)))
common_parameters(::Val{:atomic}) = ustrip.(common_parameters(Val(:atomic_unitful)))

duration(::Val{:SI_unitful}) = 18.0u"fs"
duration(::Val{:SI}) = ustrip.(u"s", duration(Val(:SI_unitful)))
duration(::Val{:atomic_unitful}) = auconvert.(duration(Val(:SI_unitful)))
duration(::Val{:atomic}) = ustrip.(duration(Val(:atomic_unitful)))
