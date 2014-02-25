module Types

export Fish, Fishers, Vars, Tau, Output

#### Define Fish type
type Fish
    xy::Array{Float64} # fish location in xy
    ci::Array{Float64} # index of nearest cluster centre to each fish
    cl::Array{Float64} # cluster xy
end

#### Define Fisher type
type Fishers
    xy::Array{Float64}
    H::Array{Int}
end

#### Define Variable type
type Vars
    Dmin::Array{Float64}
    DDx::Array{Float64}
    DDy::Array{Float64}
    ANG::Array{Float64}
    VR::Array{Float64}
    RN::Array{Float64}
    JJ::Array{Int}
    KK::Array{Int}
end

#### Define Variable type
type Tau
    n::Array{Int}
    t::Array{Int}
    s::Array{Float64}
    mu::Array{Float64}
    S::Array{Float64}
    M::Array{Float64}
    s2::Array{Float64}
    dmu::Array{Float64}
    ds2::Array{Float64}
end

#### Output variable for plotting
type Output
    fish_xy::Array{Float64}
    cons_xy::Array{Float64}
end

end
