@inline fundamental_constants(laser::AbstractLaser) = getfield(laser, :constants)
@inline fundamental_constants(laser::AbstractLaser, k) = getfield(laser.constants, k)
@inline immutable_cache(laser::AbstractLaser) = getfield(laser, :derived)
@inline immutable_cache(laser::AbstractLaser, k) = getfield(laser.derived, k)
@inline mutable_cache(laser::AbstractLaser) = getfield(laser, :cache)[]
@inline update_cache!(laser::AbstractLaser, k, v) = setfield!(laser.cache[], k, v)
@inline geometry(laser::AbstractLaser) = getfield(laser, :geometry)
@inline polarization(laser::AbstractLaser) = getfield(laser, :polarization)
@inline polarization(laser::AbstractLaser, k) = getfield(laser.polarization, k)
@inline profile(laser::AbstractLaser) = getfield(laser, :profile)

@inline function w(z, laser)
    w₀ = laser.w₀
    z_R = immutable_cache(laser, :z_R)

    w₀ * √(1 + (z/z_R)^2)
end

@inline R(z, z_R) = z + z_R^2 / z

@inline rotate_coords(R, r) = R * r
@inline rotate_coords(::UniformScaling, r) = r
@inline inv_rotate(R, r) = R \ r
@inline inv_rotate(::UniformScaling, r) = r

function Base.:(==)(a::AbstractLaser, b::AbstractLaser)
    typeof(a) == typeof(b) &&
    fundamental_constants(a) == fundamental_constants(b) &&
    immutable_cache(a) == immutable_cache(b) &&
    mutable_cache(a) == mutable_cache(b) &&
    geometry(a) == geometry(b) &&
    polarization(a) == polarization(b) &&
    profile(a) == profile(b) &&
    true
end

function Base.isequal(a::AbstractLaser, b::AbstractLaser)
    t1 = typeof(a)
    t2 = typeof(b)

    isequal(t1, t2) &&
    fundamental_constants(a) == fundamental_constants(b) &&
    immutable_cache(a) == immutable_cache(b) &&
    mutable_cache(a) == mutable_cache(b) &&
    geometry(a) == mutable_cache(b) &&
    polarization(a) == polarization(b) &&
    profile(a) == profile(b) &&
    true
end

function Base.hash(l::AbstractLaser, h::UInt)
    constants = fundamental_constants(l)
    derived = immutable_cache(l)
    cache = mutable_cache(l)
    geo = geometry(l)
    pol = polarization(l)
    prof = profile(l)

    typename = Symbol(typeof(l))
    hash(prof, hash(pol, hash(geo, hash(cache, hash(derived, hash(constants,
        hash(typename, h)))))))
end
