module Constants

export GRD_nx, GRD_dx, GRD_mx, GRD_x
export PF_N, PF_n, PF_sig
export PC_n, PC_v, PC_f, PC_j, PC_r1, PC_r2, PC_rp, PC_q, PC_h
export PS_n, PS_p

###### PARAMETERS / CONSTANTS for a seasonal run
#### Domain parameters
const GRD_nx = 100.;
const GRD_dx = 1.;
const GRD_mx = GRD_nx.*GRD_dx;
const GRD_x  = [0:GRD_dx:GRD_mx];

##### School parameters
const PS_n   = 1;    # number of schools

#### Fish parameters
const PF_n	  = 60 # number of fish per school
const PF_N 	  = PF_n * PS_n # total number of fish in the system
const PF_sig  = 3.; # distance parameter

##### Fisher parameters
const PC_n    = 6; # number of fishers
const PC_v    = .6; # speed of fishers (km per time)
const PC_h    = .5; # distance at which fishers can catch fish (km)
const PC_r1   = .1; # correlated random walk constant (steam)
const PC_r2   = .3; # correlated random walk constant (search)
const PC_rp   = 0.3; # steam / search switching probability

##### Parameters that are derived
const PS_p   = 1 ./ (4.*(GRD_mx ./ (2*PC_v))); # probability school will move
const PC_q	 = 0.05; # prob of catching fish
const PC_f   = PF_sig; # radius of fish finder (km)


end
