@doc """
    struct LaserGeometry

The geometrical aspects of the laser can be described through its orientation
with respect to a given reference system and by its propagation.
"""
LaserGeometry

struct LaserGeometry{O,R}
    orientation::O
    rotation_matrix::R
end

function LaserGeometry(orientation; isversor=false)
    v₁ = !isversor ? normalize(orientation[1]) : orientation[1]
    v₃ = !isversor ? normalize(orientation[2]) : orientation[2]
    v₂ = v₃ × v₁

    # original basis
    e₁ = SVector{3}(1,0,0)
    e₂ = SVector{3}(0,1,0)
    e₃ = SVector{3}(0,0,1)

    # original basis matrix
    E = [e₁ e₂ e₃]

    # desired direction matrix
    V = [v₁ v₂ v₃]

    # rotation matrix
    R = E * transpose(V)

    LaserGeometry(orientation, R)
end

function sym2vec(dir::Symbol)
    if dir == :x
        SVector{3}(1,0,0)
    elseif dir == :y
        SVector{3}(0,1,0)
    else
        SVector{3}(0,0,1)
    end
end

function LaserGeometry(orientation::Tuple{Symbol,Symbol})
    if orientation[1] == :x && orientation[2] == :z
        return LaserGeometry((:x, :z), I)
    end
    v₁ = sym2vec(orientation[1])
    v₃ = sym2vec(orientation[2])

    LaserGeometry((v₁, v₃), isversor=true)
end
