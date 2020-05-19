function setup_laser(laser, units; τ=nothing, kwargs...)
    _c, _q, _m_q = fundamental_constants(Val(units))
    _λ, _w₀ = common_parameters(Val(units))
    c = get(kwargs, :c, _c)
    q = get(kwargs, :q, _q)
    m_q = get(kwargs, :m_q, _m_q)
    λ = get(kwargs, :λ, _λ)
    w₀ = get(kwargs, :w₀, _w₀)

    if τ ≡ nothing
        τ = duration(Val(units))
    end
    profile = get(kwargs, :profile, GaussProfile(c=c, τ=τ))

    others = Dict{Symbol,Any}()
    excluded = [:c, :q, :m_q, :λ, :w₀, :profile]
    for (k,v) in kwargs
        if k ∉ excluded
            push!(others, k=>v)
        end
    end

    laser(c=c, q=q, m_q=m_q, λ=λ, w₀=w₀, profile=profile, others...)
end

fundamental_constants(::Val{:SI_unitful}) = c_0, -e, m_e
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
