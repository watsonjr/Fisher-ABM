#module Constants

#export PRM, set_constants

###### PARAMETERS / CONSTANTS for a seasonal run
#### Domain parameters
GRD_nx = 100; # grid size (number)
GRD_dx = 1; # grid cell size (km)
GRD_mx = GRD_nx.*GRD_dx; # grid size (km)
GRD_mx2 = GRD_mx / 2 # half grid size (used in periodic bnd)
GRD_x  = [0:GRD_dx:GRD_mx];

##### School parameters
PS_n   = 1 # number of schools
PS_p   = 0.001; # probability school will move (#do I need this)

#### Fish parameters
PF_n	 = 60 # number of fish per school
PF_sig = 5; # distance parameter (km)

##### Fisher parameters
PC_n   = 1; # number of fishers
PC_v   = 4; # max speed of fishers (km per time)
PC_h   = 1; # distance at which fishers can catch fish (km)
#const PC_r   = .25; # correlated random walk angle (drunk ballistic)
PC_rp  = 0.99; # choose random change in walk (#previously: probability of r1->r2 (r2->r1 = 1-PC_rp))
PC_f   = 5.; # radius of fish finder (x grid cells; km)
PC_q	 = 1.; # prob of catching fish
PC_lambda = 1.; #default weight in the social networkx
PC_spy = GRD_mx; #Spying distance (important to value (infinite-distance) sharing over mutual spying)


type Param
    GRD_nx::Int
    GRD_dx::Float64
    GRD_mx::Float64
    GRD_mx2::Float64
    GRD_x::Array{Int}
    PS_n::Int 
    PS_p::Float64
    PF_n::Int
    PF_sig::Float64
    PC_n::Int
    PC_v::Float64
    PC_h::Float64
    PC_rp::Float64
    PC_f::Float64
    PC_q::Float64
    PC_lambda::Float64
    PC_spy::Float64
end


PRM=Param(
GRD_nx, GRD_dx, GRD_mx, GRD_mx2, GRD_x, 
PS_n, PS_p, PF_n, PF_sig,
 PC_n, PC_v, PC_h, PC_rp, PC_f,  PC_q, PC_lambda, PC_spy)



macro set_constants(struct)
    ### This macro is black magic.
    ### Its purpose is to insert PS_n = PRM.PS_n and so on
    ### automatically for all the fields of PRM the parameter container
    ### at the beginning of every single function in the program
    #
    # We want to build up a block of expressions.
    block = Expr(:block)
    for f in [:GRD_nx,:GRD_dx,:GRD_mx,:GRD_mx2,:GRD_x,
            :PS_n, :PS_p, :PF_n, :PF_sig,
            :PC_n, :PC_v, :PC_h, :PC_rp, :PC_f,
            :PC_q, :PC_lambda, :PC_spy]
        # each expression in the block consists of
        # the fieldname = struct.fieldname
        e = :($f = $struct.$f)
        # add this new expression to our block
        push!(block.args, e)
    end
    # now escape the evaled block so that the
    # new variable declarations get declared in the surrounding scope.
    return esc(:($block))
end
#end
