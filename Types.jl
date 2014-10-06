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
    S::Array{Int} # make(1)/break(-1) friendships
    MI::Array{Int} # index of steaming or searching
    SN::Array{Float64} # social network
    V::Array{Float64} # speed
    measure::Dict{ASCIIString,Array{Float64}} # every quantity (times, frequencies) we want to measure
    #Ts::Array{Float64} # running mean time between schools for each fisher
    #Tv::Array{Float64} # running variance time between schools for each fisher
    #ts::Array{Float64} # current time between schools for each fisher
    #ns::Array{Int} # number of schools visited

end

#### Output variable for plotting
type Output
    fish_xy::Array{Float64}
    cons_xy::Array{Float64}
    schl_xy::Array{Float64}
    cons_H::Array{Float64}
end

#### Rtree for localizing nearest neighbor among fish
type Fishtree
    trees::Array{PyObject}
end

end

