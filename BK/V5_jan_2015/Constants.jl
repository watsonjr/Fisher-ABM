#module Constants

#export PRM, set_constants

###### PARAMETERS / CONSTANTS for a seasonal run
#### Domain parameters
GRD_nx = 100; # grid size (number)
GRD_dx = 1; # grid cell size (km)
GRD_mx = GRD_nx.*GRD_dx; # grid size (km)
GRD_mx2 = GRD_mx / 2 # half grid size (used in periodic bnd)

#### Fish parameters
PF_n   = 20 # number of fish per school
PF_sig = GRD_mx*0.035; # radius of school

##### Fisher parameters
PC_n   = 1; # number of fishers
PC_h   = 1; # distance at which fishers can catch fish (km)
PC_v   = PF_sig; # max speed of fishers (km per time)

PC_rp  = 0.5; # Fraction of time spent in ballistic mode

PC_f   = GRD_mx*0.0275; # radius of fish finder (x grid cells; km)
PC_q   = 1; # prob of catching fish (or rate of depletion if implicit fish)

PC_lambda = 0.; #default weight in the social network
PC_ncliq = 1; #number of cliques in the social network
PC_spy = GRD_mx; #Spying distance 

##### School parameters
PS_n   = int(round(.03/(pi*PF_sig^2/GRD_nx^2))); # number of schools
PS_p   = 0.0; # probability school will move (#do I need this)

type Param
    GRD_nx::Int
    GRD_dx::Float64
    GRD_mx::Float64
    GRD_mx2::Float64
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
    PC_ncliq::Int
    PC_spy::Float64
end


PRM=Param(
GRD_nx, GRD_dx, GRD_mx, GRD_mx2, 
PS_n, PS_p, PF_n, PF_sig,
 PC_n, PC_v, PC_h, PC_rp, PC_f,  PC_q, PC_lambda,PC_ncliq, PC_spy)

DFT_PRM=deepcopy(PRM) #Default parameters

