# shorten float strings
xprint(x) = sprint(show, x, context = :compact => true) 
# return plain name of object
wrapper_name(x) = typeof(x).name.wrapper

# general single laser verbose string
function laser_description(s)
    unit_system = fundamental_constants(s, :unit_system)
    rows = String[]
    push!(rows, "$(wrapper_name(s)) with $unit_system units")
    a_0 = xprint(s.a₀)
    push!(rows, "a₀ = "*a_0)
    omega  = xprint(immutable_cache(s, :ω)) 
    push!(rows, "ω  = "*omega)
    phi_0 = xprint(s.ϕ₀) 
    push!(rows, "ϕ₀ = "*phi_0)
    push!(rows, "and temporal $(s.profile)")
    return rows
end

# Abstract Laser
function Base.show(io::IO, ::MIME"text/plain", s::AbstractLaser)
    laser_string = join(laser_description(s),"\n")
    return print(io, laser_string)
end

# Gaussian Laser
function Base.show(io::IO, ::MIME"text/plain", s::GaussLaser)
    rows = laser_description(s)
    lastrow = rows[end]
    w_0 = xprint(s.w₀)
    rows[end] = "w₀ = "*w_0
    push!(rows, "with $(s.polarization)")
    push!(rows, lastrow)
    laser_string = join(rows, "\n")
    return print(io, laser_string)
end

# Laguerre-Gauss Laser
function Base.show(io::IO, ::MIME"text/plain", s::LaguerreGaussLaser)
    rows = laser_description(s)
    lastrow = rows[end]
    w_0 = xprint(s.w₀)
    rows[end] = "w₀ = "*w_0
    push!(rows, "p = $(s.p) and m = $(s.m)")
    push!(rows, "with $(s.polarization)")
    push!(rows, lastrow)
    laser_string = join(rows, "\n")
    return print(io, laser_string)
end

# general single laser brief string
function laser_brief(s)
    a_0 = xprint(s.a₀)
    omega = xprint(immutable_cache(s, :ω))
    phi_0 = xprint(s.ϕ₀) 
    temp = wrapper_name(s.profile)
    rows = String[]
    push!(rows, "a₀ = $a_0")
    push!(rows, "ω = $omega")
    push!(rows, "ϕ₀ = "*phi_0)
    push!(rows, "temporal $temp")
    return rows
end

# multiobject/brief print
# AbstractLaser
function Base.show(io::IO, s::AbstractLaser) 
    name = wrapper_name(s)
    rows = laser_brief(s)
    laser_string = join(rows, ", ")
    return print(io, "$name: "*laser_string)    
end 

# GaussLaser
function Base.show(io::IO, s::GaussLaser) 
    rows = laser_brief(s)
    lastrow = rows[end]
    rows[end] = "w₀ = $(xprint(s.w₀))"
    push!(rows, lastrow)
    laser_string = join(rows, ", ")
    return print(io, "Gaussian laser: "*laser_string)
end 

# LaguerreGaussLaser
function Base.show(io::IO, s::LaguerreGaussLaser) 
    rows = laser_brief(s)
    lastrow = rows[end]
    rows[end] = "w₀ = $(xprint(s.w₀))"
    push!(rows, "p = $(s.p)")
    push!(rows, "m = $(s.m)")
    push!(rows, lastrow)
    laser_string = join(rows, ", ")
    return print(io, "Laguerre-Gauss laser: "*laser_string)
end 

# polarization show
function Base.show(io::IO, p::LaserPolarization)
    @unpack ξx, ξy = p
    tangent = ξy/ξx
    if abs(ξx) ≈ abs(ξy)
        angle(tangent) ≈  π/2 && return print(io, "right-handed polarization") 
        angle(tangent) ≈ -π/2 && return print(io, "left-handed polarization") 
    else 
        angle(tangent) ≈  0.  && return print(io, "linear polarization") 
        angle(tangent) ≈  π   && return print(io, "linear polarization") 
        return print(io, "eliptical polarization")
    end
end

# temporal profile show
# AbstractTemporalProfile
function Base.show(io::IO, p::AbstractTemporalProfile)
    print(io, wrapper_name(p))
end

# GaussProfile
function Base.show(io::IO, p::GaussProfile)
    @unpack inv_τ, t₀, z₀ = p
    rows = String[]
    push!(rows, "Gaussian profile")
    z0, t0 = xprint(z₀), xprint(t₀)
    push!(rows, "centered in z₀ = $z0 and t₀ = $t0")
    tau = xprint(inv(inv_τ))
    push!(rows, "with duration of pulse (FWHM) τ = $tau")
    laser_string = join(rows, "\n")
    return print(io, laser_string)
end

# Cos²Profile
function Base.show(io::IO, p::Cos²Profile)
    @unpack τ, t₀, z₀ = p
    rows = String[]
    push!(rows, "Cos² profile")
    z0, t0 = xprint(z₀), xprint(t₀)
    push!(rows, "centered in z₀ = $z0 and t₀ = $t0")
    tau = xprint(τ)
    push!(rows, "with duration of pulse τ = $tau")
    laser_string = join(rows, "\n")
    return print(io, laser_string)
end

# QuasiRectangularProfile
function Base.show(io::IO, p::QuasiRectangularProfile)
    @unpack Δt, inv_τ, t₀, z₀ = p
    rows = String[]
    push!(rows, "Quasi rectangular profile")
    z0, t0 = xprint(z₀), xprint(t₀)
    push!(rows, "centered in z₀ = $z0 and t₀ = $t0")
    tau = xprint(inv(inv_τ))
    push!(rows, "with exponential decay width τ = $tau")
    dt = xprint(Δt)
    push!(rows, "length of the plateau Δt = $dt")
    laser_string = join(rows, "\n")
    return print(io, laser_string)
end
