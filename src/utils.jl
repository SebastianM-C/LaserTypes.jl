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

# The cache field is excluded from equality and hashing: it holds transient
# task-local scratch state, so including it would make the results depend on
# the evaluation history of the current task.
function Base.:(==)(a::AbstractLaser, b::AbstractLaser)
    typeof(a) == typeof(b) || return false
    for f in fieldnames(typeof(a))
        f === :cache && continue
        getfield(a, f) == getfield(b, f) || return false
    end
    return true
end

function Base.isequal(a::AbstractLaser, b::AbstractLaser)
    typeof(a) == typeof(b) || return false
    for f in fieldnames(typeof(a))
        f === :cache && continue
        isequal(getfield(a, f), getfield(b, f)) || return false
    end
    return true
end

function Base.hash(l::AbstractLaser, h::UInt)
    h = hash(Symbol(typeof(l)), h)
    for f in fieldnames(typeof(l))
        f === :cache && continue
        h = hash(getfield(l, f), h)
    end
    return h
end
