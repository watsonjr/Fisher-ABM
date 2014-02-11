
#### Preamble
using NPZ, Distance, Devectorize
include("sub_functions.jl");
include("sub_dynamics.jl");


#### Run a season
include("make_params_s.jl");
P_C_Sn = ones(P_C_n,P_C_n);
fish_xy,cons_xy,cons_H = make_season_s(Fish_xy,Cons_xy,Cons_H,P_C_Sn);
P_C_Sn = eye(P_C_n);
ffish_xy,ccons_xy,ccons_H = make_season_s(Fish_xy,Cons_xy,Cons_H,P_C_Sn);

#### Test different social networks at equilibrium
include("make_params_e.jl");
P_C_Sn = ones(P_C_n,P_C_n);
Tau_mu1 = make_equilibrium(Fish_xy,Cons_xy,Cons_H,P_C_Sn);
P_C_Sn = eye(P_C_n);
Tau_mu2 = make_equilibrium(Fish_xy,Cons_xy,Cons_H,P_C_Sn);


#### Save
npzwrite("./Data/Data_fish.npy", Fish_xy)
npzwrite("./Data/Data_fisher_xy.npy", Cons_xy)
npzwrite("./Data/Data_fisher_H.npy", Cons_H)
npzwrite("./Data/Data_fclust.npy", Fclust_xy)

#### Useful Julia code
###### Grid fish density
#ii = int(round(Fish_xy[:,:,1]));
#jj,mm,kk = hist2d(ii,Grd_x[1:2:end],Grd_y[1:2:end]);

