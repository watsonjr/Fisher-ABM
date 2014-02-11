###### Parameters for fish agent-based model

#### Domain parameters
Grd_nx = 100.;
Grd_ny = 100.;
Grd_dx = 1.;
Grd_dy = 1.;
Grd_mx = Grd_nx.*Grd_dx;
Grd_my = Grd_ny.*Grd_dy;
Grd_x  = [0:Grd_dx:Grd_mx];
Grd_y  = [0:Grd_dy:Grd_my];

#### Model parameters
P_F_n    = 200; # number of fish
P_F_dx   = 2; # distance fish are placed from cluster centre
P_F_m    = 0.2; # max fraction of fish that relocate each time step

P_C_n    = 5; # number of fishers
P_C_v    = 0.5; # speed of fishers
P_C_dx   = 5; # radius of fish finder
P_C_Sn   = eye(P_C_n); # fisher social network

P_Cl_n   = 5;    # number of cluster centres
P_Cl_ang = 30;   # cluster angle change per unit time (i.e. +- Fc_ang)
P_Cl_al  = 1.5;  # levy flight constant (determins max dist moved per unit time; ref 1)
P_Cl_dx  = .01;  # min distance that clusters can move per unit time

P_Tend  = 1000; # number of timesteps
P_ddx   = 10.; # P_ddx/2 = the dx and dy of fish finder
P_catch = 4; # distance at which fishers can catch fish
P_prob  = 0.5; # prob## Stuff I use
P_DXY   = zeros(P_C_n,2); # direction vector used by fishers
P_vr    = zeros(P_C_n);  # speed modulator (proportional to distance to fish)
#theta_f = zeros(F_n,1);

##### Initialize clusters
Fclust_xy = zeros(P_Cl_n,2,P_Tend); # cluster location through time
Fclust_xy[:,:,1]=rand(Float64,P_Cl_n,2).*repmat([Grd_nx*Grd_dx Grd_ny*Grd_dy],P_Cl_n,1);

##### Initial fish
Fish_xy = zeros(P_F_n,2,P_Tend);
for i = 1:P_F_n # place each fish
    j = ceil(rand(1)*P_Cl_n); # index of which cluster center
    Fish_xy[i,:,1] = mod(Fclust_xy[j,:,1] + randn(1,2).*P_F_dx, [Grd_mx Grd_my]);
end

##### Initialize fishers
Cons_H     = zeros(P_C_n,P_Tend); # catch at time t
Cons_xy    = zeros(P_C_n,2,P_Tend); # Fisher locations through time
Cons_xy[:,:,1] = rand(Float64,P_C_n,2).*repmat([Grd_nx*Grd_dx Grd_ny*Grd_dy],P_C_n,1);



#### REFS
# 1 - http://math.stackexchange.com/questions/52869/numerical-approximation-of-levy-flight/52870#52870
