using LaserTypes
using LaserTypes: pochhammer
using CSV
using Test

ω = 0.057
c = 137.036
λ = 2π*c/ω
w₀ = 75*λ
s = setup_laser(LaguerreGaussLaser, :atomic)
ξx_dict = Dict("1/Sqrt[2]"=>1/√2, "1"=>1.)
ξy_dict = Dict("1/Sqrt[2]"=>im/√2, "0"=>0., "(-I)/Sqrt[2]"=>-im/√2)

@testset "Comparison to Symbolic Values" begin
    for row in CSV.File("field.csv", header = false)
        if row[begin] == "p"
            p = convert(Int64, row[2])
            m = convert(Int64, row[4])
            ξx = ξx_dict[row[6]]
            ξy = ξy_dict[row[8]]
            s = setup_laser(LaguerreGaussLaser, :atomic; λ, w₀, ξx, ξy, p, m, profile = ConstantProfile())
        else
            x = tryparse(Float64, row[1])
            y = tryparse(Float64, row[2])
            z = tryparse(Float64, row[3])
            EX = tryparse(ComplexF64, replace(replace(row[4],"*I"=>"im"),"*^"=>"e"))
            EY = tryparse(ComplexF64, replace(replace(row[5],"*I"=>"im"),"*^"=>"e"))
            EZ = tryparse(ComplexF64, replace(replace(row[6],"*I"=>"im"),"*^"=>"e"))
            BX = tryparse(ComplexF64, replace(replace(row[7],"*I"=>"im"),"*^"=>"e"))
            BY = tryparse(ComplexF64, replace(replace(row[8],"*I"=>"im"),"*^"=>"e"))
            BZ = tryparse(ComplexF64, replace(replace(row[9],"*I"=>"im"),"*^"=>"e"))
            juliaE = E([x,y,z], s)
            juliaB = B([x,y,z], s)
            @test juliaE[1] ≈ factorial(s.p)/pochhammer(abs(s.m)+1,s.p)*s.E₀*s.Nₚₘ*EX
            @test juliaE[2] ≈ factorial(s.p)/pochhammer(abs(s.m)+1,s.p)*s.E₀*s.Nₚₘ*EY
            @test juliaE[3] ≈ factorial(s.p)/pochhammer(abs(s.m)+1,s.p)*s.E₀*s.Nₚₘ*EZ
            @test juliaB[1] ≈ factorial(s.p)/pochhammer(abs(s.m)+1,s.p)*s.E₀*s.Nₚₘ*BX
            @test juliaB[2] ≈ factorial(s.p)/pochhammer(abs(s.m)+1,s.p)*s.E₀*s.Nₚₘ*BY
            @test juliaB[3] ≈ factorial(s.p)/pochhammer(abs(s.m)+1,s.p)*s.E₀*s.Nₚₘ*BZ
        end
    end
end
