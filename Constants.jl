module Constants

export GRD_nx, GRD_dx, GRD_mx, GRD_x
export PF_N, PF_n, PF_sig
export PC_n, PC_v, PC_f, PC_j, PC_r1, PC_r2, PC_rp, PC_q, PC_h
export PS_n, PS_p

###### PARAMETERS / CONSTANTS for a seasonal run
#### Domain parameters
const GRD_nx = 100;
const GRD_dx = 1;
const GRD_mx = GRD_nx.*GRD_dx;
const GRD_x  = [0:GRD_dx:GRD_mx];

##### School parameters
const PS_n   = 2 # number of schools

#### Fish parameters
const PF_n	 = 60 # number of fish per school
const PF_N 	 = PF_n * PS_n # total number of fish in the system
const PF_sig = 3.; # distance parameter (km)

##### Fisher parameters
const PC_n   = 6; # number of fishers
const PC_v   = 1.2; # speed of fishers (km per time)
const PC_h   = .5; # distance at which fishers can catch fish (km)
const PC_r1  = .05; # correlated random walk constant (steam)
const PC_r2  = .4; # correlated random walk constant (search)
const PC_rp  = 0.15; # steam / search switching probability
const PC_f   = 5.; # radius of fish finder (x grid cells; km)

##### Parameters that are derived
const PS_p   = 1 ./ (4.*(GRD_mx ./ (2*PC_v))); # probability school will move
const PC_q	 = .5; # prob of catching fish


end
