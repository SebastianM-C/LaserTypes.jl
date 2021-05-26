using LaserTypes
using LaserTypes: pochhammer, immutable_cache
using CSV
using UnPack
using Test

ω = 0.057
c = 137.03599908330935
λ = 2π*c/ω
w₀ = 75*λ

normalization(m, p, E₀, Nₚₘ) = factorial(p)/pochhammer(abs(m)+1,p)*E₀*Nₚₘ

csv_folder = joinpath(@__DIR__, "fields_csv/")

for (roots, dirs, files) in walkdir(csv_folder)
    @testset "$file" for file in files
        i, p, m = tryparse.(Int, (split(file,"_"))[2:4])
        ξx = cos(i*π/4) + 0im; ξy = im*sin(i*π/4);
        s = setup_laser(LaguerreGaussLaser, :atomic; λ, w₀, p, m, ξx, ξy, profile = ConstantProfile())
        @unpack E₀, Nₚₘ = immutable_cache(s)

        for row in CSV.File(joinpath(roots, file), header=false)
            x = row[1]
            y = row[2]
            z = row[3]

            EX = tryparse(ComplexF64, replace(replace(row[4],"*I"=>"im"),"*^"=>"e"))
            EY = tryparse(ComplexF64, replace(replace(row[5],"*I"=>"im"),"*^"=>"e"))
            EZ = tryparse(ComplexF64, replace(replace(row[6],"*I"=>"im"),"*^"=>"e"))

            BX = tryparse(ComplexF64, replace(replace(row[7],"*I"=>"im"),"*^"=>"e"))
            BY = tryparse(ComplexF64, replace(replace(row[8],"*I"=>"im"),"*^"=>"e"))
            BZ = tryparse(ComplexF64, replace(replace(row[9],"*I"=>"im"),"*^"=>"e"))

            juliaE = E([x,y,z], s)
            juliaB = B([x,y,z], s)

            @test juliaE[1] ≈ normalization(s.m, s.p, E₀, Nₚₘ)*EX
            @test juliaE[2] ≈ normalization(s.m, s.p, E₀, Nₚₘ)*EY
            @test juliaE[3] ≈ normalization(s.m, s.p, E₀, Nₚₘ)*EZ

            @test juliaB[1] ≈ normalization(s.m, s.p, E₀, Nₚₘ)*BX
            @test juliaB[2] ≈ normalization(s.m, s.p, E₀, Nₚₘ)*BY
            @test juliaB[3] ≈ normalization(s.m, s.p, E₀, Nₚₘ)*BZ
        end
    end
end
