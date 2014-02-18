
############ typeS
### Define parameter type
immutable T_Params
    G_nx::Float64
    G_ny::Float64
    G_dx::Float64
    G_dy::Float64
    G_mx::Float64
    G_my::Float64
    G_x::Array{Float64}
    G_y::Array{Float64}
    F_n::Int
    F_dx::Float64
    F_m::Float64
    C_n::Int
    C_v::Float64
    C_ff::Float64
    C_c::Float64
    C_ang::Float64
    C_pr::Float64
    Cl_n::Int
    Cl_id::Array{Int}
    Cl_ang::Float64
    Cl_al::Float64
end

#### Define Fish type
type T_Fish
    xy::Array{Float64} # fish location in xy
    ci::Array{Float64} # index of nearest cluster centre to each fish
    cl::Array{Float64} # cluster xy
end

#### Define Fisher type
type T_Fishers
    xy::Array{Float64}
    H::Array{Int}
end

#### Define Variable type
type T_Vars
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
type T_Tau
    n::Array{Int}
    t::Array{Int}
    s::Array{Float64}
    mu::Array{Float64}
    dmu::Array{Float64}
end

