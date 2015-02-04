module Types

using PyCall

export School, Fish, Fishers, Fishtree, Output

#### Define Fish type
type Fish
    fx::Array{Float64} # fish location in xy
    fs::Array{Int} # index of school fish is associated with
end

type School
    x::Array{Float64} # school xy
    fish::Array{Int} #index of the fish in the school
    pop::Array{Float64} #population of each school (for implicit fish case)
end

#### Define Fisher type
type Fishers
    x::Array{Float64} # location
    Ni::Array{Int} # index of nearest fish
    target::Array{Int} # index of target fish (or school)
    Dmin::Array{Float64} # distance to target fish
    DXY::Array{Float64} # direction unit vector
    H::Array{Float64} # harvest count
    MI::Array{Int} # index of steaming or searching
    SN::Array{Float64} # social network
    V::Array{Float64} # speed
    #HM::Array{Float64} # homophili network (number of times worked with someone)
    measure::Dict{ASCIIString,Array{Float64}} # Things to measure
end

#### Output variable for plotting
type Output
    fish_xy::Array{Float64}
    cons_xy::Array{Float64}
    schl_xy::Array{Float64}
    schl_pop::Array{Float64}
    cons_H::Array{Float64}
end

#### Rtree for localizing nearest neighbor among fish
type Fishtree
    trees::Array{PyObject}
end

end

