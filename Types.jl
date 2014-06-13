module Types

export Fish, Fishers, Output

#### Define Fish type
type Fish
    fx::Array{Float64} # fish location in xy
    fs::Array{Int} # index of school fish is associated with
    sx::Array{Float64} # school xy
end

#### Define Fisher type
type Fishers
     x::Array{Float64} # location
    H::Array{Float64} # harvest count
    S::Array{Int} # make(1)/break(-1) friendships
    MI::Array{Int} # index of steaming or searching
    CN::Array{Int} # contact network
    Dmin::Array{Float64} # distance to nearest fish
    DXY::Array{Float64} # direction vector (unit)
    VR::Array{Float64} # speed
    JJ::Array{Int} # index of nearest fish
    KK::Array{Int} # 1/0 harvest index
    cs::Array{Float64} # cumulative harvest
    Dist::Array{Float64} # cumulative distance traveled 
end

#### Output variable for plotting
type Output
    fish_xy::Array{Float64}
    cons_xy::Array{Float64}
    schl_xy::Array{Float64}
    cons_H::Array{Float64}
end

end
