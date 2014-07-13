module Constants

export GRD_nx, GRD_dx, GRD_mx, GRD_x, GRD_mx2
export PF_N, PF_n, PF_sig
export PC_n, PC_v, PC_f, PC_j, PC_r1, PC_r2, PC_rp, PC_q, PC_h, PC_vmn
export PS_n, PS_p

###### PARAMETERS / CONSTANTS for a seasonal run
#### Domain parameters
const GRD_nx = 100;
const GRD_dx = 1;
const GRD_mx = GRD_nx.*GRD_dx;
const GRD_mx2 = GRD_mx / 2
const GRD_x  = [0:GRD_dx:GRD_mx];

##### School parameters
const PS_n   = 1 # number of schools
const PS_p   = 0.001; # probability school will move (#do I need this)

#### Fish parameters
const PF_n	 = 60 # number of fish per school
const PF_N 	 = PF_n * PS_n # total number of fish in the system
const PF_sig = 2.; # distance parameter (km)

##### Fisher parameters
const PC_n   = 1; # number of fishers
const PC_v   = 2.6; # max speed of fishers (km per time)
const PC_vmn = PC_v/5; # min speed of fishers (km per time)
const PC_h   = 1; # distance at which fishers can catch fish (km)
const PC_r1  = .005; # correlated random walk constant (tumble)
const PC_r2  = 1.; # correlated random walk constant (ballistic)
#const PC_rp  = 0.01; # probability of r1->r2 (r2->r1 = 1-PC_rp)
const PC_f   = 5.; # radius of fish finder (x grid cells; km)
const PC_q	 = 1.; # prob of catching fish


end
