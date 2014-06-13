module Constants

export GRD_nx, GRD_ny, GRD_dx, GRD_dy, GRD_mx, GRD_my, GRD_x, GRD_y
export PF_n, PF_m, PF_g, PF_dx, PF_DX
export PC_n, PC_v, PC_ff, PC_c, PC_rn, PC_q, PC_Q
export PCL_n, PCL_p

###### PARAMETERS / CONSTANTS for a seasonal run
#### Domain parameters
const GRD_nx = 100.;
const GRD_ny = 100.;
const GRD_dx = 1.;
const GRD_dy = 1.;
const GRD_mx = GRD_nx.*GRD_dx;
const GRD_my = GRD_ny.*GRD_dy;
const GRD_x  = [0:GRD_dx:GRD_mx];
const GRD_y  = [0:GRD_dy:GRD_my];

##### School parameters
const PCL_n   = 2;    # number of cluster centres
const PCL_p   = 0.02; # probability that cluster center will move

#### Fish parameters
const PF_N    = 100; # expected total number of fish
const PF_n	  = PF_N / PS_n; # expected number of fish per school
const PF_m    = 0.0; # max fraction of fish that relocate each time step
const PF_g    = 0.5; # probability of fish regrowth
const PF_dx   = 0; # distance parameter

##### Fisher parameters
const PC_n    = 2; # number of fishers
const PC_v    = .3; # speed of fishers (km per time)
const PC_ff   = 8; # radius of fish finder (km)
const PC_c    = 1; # distance at which fishers can catch fish (km)
const PC_rn   = .1; # correlated random walk randomness
const PC_q	  = 0; # prob of catching fish once encountered

##### Parameters to sweep through
const PF_DX = linspace(0.01,1.5,10); # fish distribution variance
const PC_Q  = linspace(0.001,0.1,10);   # probability of catching fish

end
