function get_units(unit_str)
    # workaround stuff like kg.m/s
    unit_str = replace(unit_str, "."=>"*")
    # workaround stuff like 1/m^3
    if !occursin(r"1\/([a-z]*\^[0-9]?)", unit_str)
        units = @eval @u_str $unit_str
    else
        u = replace(unit_str, r"1\/([a-z]*)\^([0-9]?)" => s"\g<1>^-\g<2>")
        units = @eval @u_str $u
    end
end
