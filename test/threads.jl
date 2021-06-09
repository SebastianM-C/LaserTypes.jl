using LaserTypes
using LaserTypes: pochhammer, immutable_cache
using CSV
using UnPack
using Test
using Base.Threads

ω = 0.057
c = 137.03599908330935
λ = 2π*c/ω
w₀ = 75*λ

normalization(m, p, E₀, Nₚₘ) = factorial(p)/pochhammer(abs(m)+1,p)*E₀*Nₚₘ

csv_folder = joinpath(@__DIR__, "Laguerre-Gauss/fields_csv/")
ntrds = nthreads()

@testset "Thread safety" begin
    for (roots, dirs, files) in walkdir(csv_folder)
        for file in files
            i, p, m = tryparse.(Int, (split(file,"_"))[2:4])
            ξx = cos(i*π/4) + 0im; ξy = im*sin(i*π/4);
            s = setup_laser(LaguerreGaussLaser, :atomic; λ, w₀, p, m, ξx, ξy, profile = ConstantProfile())
            @unpack E₀, Nₚₘ = immutable_cache(s)

            current_file = CSV.File(joinpath(roots, file), header = false)
            num_of_rows = length(current_file)
            test_bool = Atomic{Bool}(true)

            
            current_row = zeros(Int, ntrds)
            @threads for id in 1:ntrds
                while current_row[id] + id ≤ num_of_rows
                    x = current_file[current_row[id] + id][1]
                    y = current_file[current_row[id] + id][2]
                    z = current_file[current_row[id] + id][3]

                    EX = tryparse(ComplexF64, replace(replace(current_file[current_row[id] + id][4],"*I"=>"im"),"*^"=>"e"))
                    EY = tryparse(ComplexF64, replace(replace(current_file[current_row[id] + id][5],"*I"=>"im"),"*^"=>"e"))
                    EZ = tryparse(ComplexF64, replace(replace(current_file[current_row[id] + id][6],"*I"=>"im"),"*^"=>"e"))

                    BX = tryparse(ComplexF64, replace(replace(current_file[current_row[id] + id][7],"*I"=>"im"),"*^"=>"e"))
                    BY = tryparse(ComplexF64, replace(replace(current_file[current_row[id] + id][8],"*I"=>"im"),"*^"=>"e"))
                    BZ = tryparse(ComplexF64, replace(replace(current_file[current_row[id] + id][9],"*I"=>"im"),"*^"=>"e"))

                    juliaE = E([x,y,z], s)
                    juliaB = B([x,y,z], s)
                
                    atomic_and!(test_bool, juliaE[1] ≈ normalization(s.m, s.p, E₀, Nₚₘ)*EX)
                    atomic_and!(test_bool, juliaE[2] ≈ normalization(s.m, s.p, E₀, Nₚₘ)*EY)
                    atomic_and!(test_bool, juliaE[3] ≈ normalization(s.m, s.p, E₀, Nₚₘ)*EZ)

                    atomic_and!(test_bool, juliaB[1] ≈ normalization(s.m, s.p, E₀, Nₚₘ)*BX)
                    atomic_and!(test_bool, juliaB[2] ≈ normalization(s.m, s.p, E₀, Nₚₘ)*BY)
                    atomic_and!(test_bool, juliaB[3] ≈ normalization(s.m, s.p, E₀, Nₚₘ)*BZ)
                    current_row[id] += ntrds
                end
            end
            @test test_bool.value
        end
    end
end
