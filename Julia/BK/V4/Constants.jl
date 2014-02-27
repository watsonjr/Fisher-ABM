module Constants

export GRD_nx, GRD_ny, GRD_dx, GRD_dy, GRD_mx, GRD_my, GRD_x, GRD_y
export PF_n, PF_dx, PF_m
export PC_n, PC_v, PC_ff, PC_c, PC_ang, PC_pr
export PCL_n, PCL_id, PCL_ang, PCL_al
export PTend

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

#### Fish parameters
const PF_n    = 300; # number of fish
const PF_dx   = 1; # distance fish are placed from cluster centre
const PF_m    = 0.0; # max fraction of fish that relocate each time step

##### Fisher parameters
const PC_n    = 8; # number of fishers
const PC_v    = .1; # speed of fishers
const PC_ff   = 10; # radius of fish finder
const PC_c    = .5; # distance at which fishers can catch fish
const PC_ang  = 1; # correlated random walk angle
const PC_pr   = 0.5; # probability of catching fish once encountered

##### Cluster parameters
const PCL_n   = 5;    # number of cluster centres
const PCL_id  = [1:5]; # cluster ids
const PCL_ang = 30;   # cluster angle change per unit time (i.e. +- Fc_ang)
const PCL_al  = 2.5;  # levy flight constant (ref 1)

##### Integration parameters
const PTend  = 10000; # number of timesteps

end
