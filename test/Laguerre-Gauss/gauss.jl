using LaserTypes
using Test

ω = 0.057
T₀ = 2π/ω
c = 137.036
λ = c*T₀
w₀ = 75 * λ
a₀ = 2.
LG = setup_laser(LaguerreGaussLaser, :atomic, m = 0, p = 0, λ = λ, a₀ = a₀, w₀ = w₀,
    profile=ConstantProfile())
G = setup_laser(GaussLaser, :atomic, λ = λ, a₀ = a₀, w₀ = w₀,
    profile=ConstantProfile())
for i in 1:100
    x, y, z = w₀ * (0.5 .- rand(3))
    t = T₀ * rand()
    @test Fμν([t*c, x, y, z], G) ≈ Fμν([t*c, x, y, z], LG) rtol = 1e-8
end
