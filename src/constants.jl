struct FundamentalConstants{C,Q,M,E,Mu}
    c::C
    q::Q
    mₑ::M
    ε₀::E
    μ₀::Mu
end

fundamental_constants(::Val{:SI_unitful}) = c_0, -e, m_e, ε_0, μ_0
fundamental_constants(::Val{:SI}) = ustrip.(fundamental_constants(Val(:SI_unitful)))
fundamental_constants(::Val{:atomic_unitful}) = auconvert.(fundamental_constants(Val(:SI_unitful)))
fundamental_constants(::Val{:atomic}) = ustrip.(fundamental_constants(Val(:atomic_unitful)))

function FundamentalConstants(units::Symbol)
    constants = fundamental_constants(Val(units))
    FundamentalConstants(constants...)
end
