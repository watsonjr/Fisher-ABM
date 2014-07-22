module Constants

export GRD_nx, GRD_dx, GRD_mx, GRD_x, GRD_mx2
export PF_N, PF_n, PF_sig
export PC_n, PC_v, PC_f, PC_j, PC_r, PC_rp, PC_q, PC_h
export PS_n, PS_p

###### PARAMETERS / CONSTANTS for a seasonal run
#### Domain parameters
const GRD_nx = 100; # grid size (number)
const GRD_dx = 1; # grid cell size (km)
const GRD_mx = GRD_nx.*GRD_dx; # grid size (km)
const GRD_mx2 = GRD_mx / 2 # half grid size (used in periodic bnd)
const GRD_x  = [0:GRD_dx:GRD_mx];

##### School parameters
const PS_n   = 1 # number of schools
const PS_p   = 0.01; # probability school will move (#do I need this)

#### Fish parameters
const PF_n	 = 60 # number of fish per school
const PF_sig = 2.; # distance parameter (km)

##### Fisher parameters
const PC_n   = 2; # number of fishers
const PC_v   = 2.6; # max speed of fishers (km per time)
const PC_h   = 1; # distance at which fishers can catch fish (km)
#const PC_r   = .25; # correlated random walk angle (drunk ballistic)
#const PC_rp  = 0.01; # probability of r1->r2 (r2->r1 = 1-PC_rp)
const PC_f   = 5.; # radius of fish finder (x grid cells; km)
const PC_q	 = 1.; # prob of catching fish


end
