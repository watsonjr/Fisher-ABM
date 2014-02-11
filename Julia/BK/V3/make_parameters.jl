###### PARAMETERS / CONSTANTS
#### Domain parameters
const Grd_nx = 100.;
const Grd_ny = 100.;
const Grd_dx = 1.;
const Grd_dy = 1.;
const Grd_mx = Grd_nx.*Grd_dx;
const Grd_my = Grd_ny.*Grd_dy;
const Grd_x  = [0:Grd_dx:Grd_mx];
const Grd_y  = [0:Grd_dy:Grd_my];
global Grd_nx,Grd_ny,Grd_dx,Grd_dy,Grd_mx,Grd_my,Grd_x,Grd,y

#### Fish parameters
const P_F_n    = 300; # number of fish
const P_F_dx   = 1; # distance fish are placed from cluster centre
const P_F_m    = 0.01; # max fraction of fish that relocate each time step
global P_F_n,P_F_dx,P_F_m

##### Fisher parameters
const P_C_n    = 8; # number of fishers
const P_C_v    = 1.5; # speed of fishers
const P_C_ff   = 10; # radius of fish finder
const P_C_c    = 3; # distance at which fishers can catch fish
const P_C_ang  = 1; # correlated random walk angle
const P_C_Sn   = eye(P_C_n); # fisher social network
const P_C_pr   = 0.5; # probability of catching fish once encountered
#const P_C_Sn   = ones(P_C_n,P_C_n);
global P_C_n,P_C_v,P_C_ff,P_C_c,P_C_Sn

##### Cluster parameters
const P_Cl_n   = 5;    # number of cluster centres
const P_Cl_id  = [1:5]; # cluster ids
const P_Cl_ang = 30;   # cluster angle change per unit time (i.e. +- Fc_ang)
const P_Cl_al  = 3.5;  # levy flight constant (ref 1)
global P_Cl_n,P_Cl_ang,P_Cl_al,P_Cl_dx

##### Integration parameters
const P_Tend  = 1000; # number of timesteps
global P_Tend


######## VARIABLES
Dmin  = zeros(Float64,P_C_n); # distance to nearest fish
DDx   = zeros(Float64,P_C_n); # x vector compontent to nearest fish
DDy   = zeros(Float64,P_C_n); # y vector compontent to nearest fish
ANG   = zeros(Float64,P_C_n); # fisher direction angle
VR    = zeros(Float64,P_C_n);  # speed modulator (proportional to distance to fish)
RN    = zeros(Float64,1,P_C_n); # random communication (in fnc_information)
JJ    = zeros(Int,P_C_n); # index of nearest fish (0 if nothing near)
KK    = zeros(Int,P_C_n); # index of catch [0,1]

##### Initialize clusters
Fclust_xy = zeros(Float64,P_Cl_n,2,P_Tend); # cluster location through time
Fclust_xy[:,:,1]=rand(Float64,P_Cl_n,2).*repmat([Grd_nx*Grd_dx Grd_ny*Grd_dy],P_Cl_n,1);

##### Initial fish
Fish_xy = zeros(Float64,P_F_n,2,P_Tend); # location of fish
Fish_cl = zeros(Float64,P_F_n); # index of cluster centre fishers are nearest
for i = 1:P_F_n # place each fish
    j = ceil(rand(1)*P_Cl_n); # index of which cluster center
    Fish_cl[i] = j[1];
    Fish_xy[i,:,1] = mod(Fclust_xy[j,:,1] + randn(1,2).*P_F_dx, [Grd_mx Grd_my]);
end

##### Initialize fishers
Cons_H     = zeros(Int,P_C_n,P_Tend); # catch at time t
Cons_xy    = zeros(Float64,P_C_n,2,P_Tend); # Fisher locations through time
Cons_xy[:,:,1] = rand(Float64,P_C_n,2).*repmat([Grd_nx*Grd_dx Grd_ny*Grd_dy],P_C_n,1);


#### REFS
# 1 - http://math.stackexchange.com/questions/52869/numerical-approximation-of-levy-flight/52870#52870
