struct LaserGeometry{D,R}
    oscillation_dir::D
    propagation_dir::D
    rotation_matrix::R
end

function LaserGeometry(oscillation_dir, propagation_dir; isversor=false)
    v₁ = !isversor ? normalize(oscillation_dir) : oscillation_dir
    v₃ = !isversor ? normalize(propagation_dir) : propagation_dir
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

    LaserGeometry(oscillation_dir, propagation_dir, R)
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

function LaserGeometry(oscillation_dir::Symbol, propagation_dir::Symbol)
    if oscillation_dir == :x && propagation_dir == :z
        return LaserGeometry(:x, :z, I)
    end
    v₁ = sym2vec(oscillation_dir)
    v₃ = sym2vec(propagation_dir)

    LaserGeometry(v₁, v₃, isversor=true)
end
