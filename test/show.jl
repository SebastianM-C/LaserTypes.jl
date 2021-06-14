using LaserTypes, Test

GCL  = setup_laser(GaussLaser, :atomic; profile = ConstantProfile())
GGL  = setup_laser(GaussLaser, :atomic; profile = GaussProfile(z₀ = 0.0, t₀ = 0.0, τ = 744.145))
G2L  = setup_laser(GaussLaser, :atomic; profile = Cos²Profile(z₀ = 0.0, t₀ = 0.0, τ = 744.145))
GQL  = setup_laser(GaussLaser, :atomic; profile = QuasiRectangularProfile(z₀ = 0.0, t₀ = 0.0, τ = 744.145, Δt = 7441.45))

LGCL = setup_laser(LaguerreGaussLaser, :atomic; profile = ConstantProfile(), p = 1, m = 1)
LGCD = setup_laser(LaguerreGaussLaser, :atomic; profile = ConstantProfile(), p = 1, m = 1, ξx = 1/√2, ξy =  1im/√2)
LGCS = setup_laser(LaguerreGaussLaser, :atomic; profile = ConstantProfile(), p = 1, m = 1, ξx = 1/√2, ξy = -1im/√2)

vec = [GCL, LGCL]

# DO NOT DISTURB THE STRINGS UNLESS YOU KNOW WHAT YOU ARE DOING OR YOU HAVE TIME TO TEST THINGS OUT!

GCL_str = 
"""GaussLaser with atomic units
a₀ = 1.0
ω  = 0.0569542
ϕ₀ = 0.0
w₀ = 1.09604e6
with linear polarization
and temporal ConstantProfile"""

GGL_str = 
"""GaussLaser with atomic units
a₀ = 1.0
ω  = 0.0569542
ϕ₀ = 0.0
w₀ = 1.09604e6
with linear polarization
and temporal Gaussian profile
centered in z₀ = 0.0 and t₀ = 0.0
with duration of pulse (FWHM) τ = 744.145"""

G2L_str =
"""GaussLaser with atomic units
a₀ = 1.0
ω  = 0.0569542
ϕ₀ = 0.0
w₀ = 1.09604e6
with linear polarization
and temporal Cos² profile
centered in z₀ = 0.0 and t₀ = 0.0
with duration of pulse τ = 744.145"""

GQL_str =
"""GaussLaser with atomic units
a₀ = 1.0
ω  = 0.0569542
ϕ₀ = 0.0
w₀ = 1.09604e6
with linear polarization
and temporal Quasi rectangular profile
centered in z₀ = 0.0 and t₀ = 0.0
with exponential decay width τ = 744.145
length of the plateau Δt = 7441.45"""

LGCL_str = 
"""LaguerreGaussLaser with atomic units
a₀ = 1.0
ω  = 0.0569542
ϕ₀ = 0.0
w₀ = 1.09604e6
p = 1 and m = 1
with linear polarization
and temporal ConstantProfile"""

LGCD_str = 
"""LaguerreGaussLaser with atomic units
a₀ = 1.0
ω  = 0.0569542
ϕ₀ = 0.0
w₀ = 1.09604e6
p = 1 and m = 1
with right-handed polarization
and temporal ConstantProfile"""

LGCS_str = 
"""LaguerreGaussLaser with atomic units
a₀ = 1.0
ω  = 0.0569542
ϕ₀ = 0.0
w₀ = 1.09604e6
p = 1 and m = 1
with left-handed polarization
and temporal ConstantProfile"""

vec_str = 
"""
2-element Vector{LaserTypes.AbstractLaser}:
 Gaussian laser: a₀ = 1.0, ω = 0.0569542, ϕ₀ = 0.0, w₀ = 1.09604e6, temporal ConstantProfile
 Laguerre-Gauss laser: a₀ = 1.0, ω = 0.0569542, ϕ₀ = 0.0, w₀ = 1.09604e6, p = 1, m = 1, temporal ConstantProfile"""

@testset "Printing/Show" begin
    io = IOBuffer()

    # Gauss Tests
    # Constant Profile
    show(io, MIME"text/plain"(), GCL)
    @test String(take!(io)) == GCL_str
    
    # Gaussian Profile
    show(io, MIME"text/plain"(), GGL)
    @test String(take!(io)) == GGL_str
    
    # Cos² Profile
    show(io, MIME"text/plain"(), G2L)
    @test String(take!(io)) == G2L_str
    
    # Quasi Rectangular Profile
    show(io, MIME"text/plain"(), GQL)
    @test String(take!(io)) == GQL_str
    
    # Laguerre Gauss Laser
    # Linearly Polarized
    show(io, MIME"text/plain"(), LGCL)
    @test String(take!(io)) == LGCL_str
    
    # Right Circularly Polarized
    show(io, MIME"text/plain"(), LGCD)
    @test String(take!(io)) == LGCD_str
    
    # Left Circularly Polarized
    show(io, MIME"text/plain"(), LGCS)
    @test String(take!(io)) == LGCS_str

    # Brief Printing Test
    # Vector Test
    show(io, MIME"text/plain"(), vec)
    @test String(take!(io)) == vec_str
end
