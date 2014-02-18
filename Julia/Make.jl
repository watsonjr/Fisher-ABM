
#### Preamble
using NPZ, Distance
include("sub_dynamics.jl");
include("sub_functions.jl");
include("sub_params.jl");

#### Initialize
include("make_types.jl");
P,Fish,Con,Var,Tau = initialize_e(); # equilibrium run

SN1 = eye(P.C_n);
SN2 = ones(P.C_n,P.C_n);

#### Run model
f = T_Fish(Fish.xy,Fish.ci,Fish.cl);
c = T_Fishers(Con.xy,Con.H);
t = T_Tau(Tau.n,Tau.t,Tau.s,Tau.mu,Tau.dmu);
v = T_Vars(Var.Dmin,Var.DDx,Var.DDy,Var.ANG,Var.VR,Var.RN,Var.JJ,Var.KK);
make_equilibrium(f,c,t,v,SN1);
out1 = t.mu;

f = Fish; c = Con; v = Var; t = Tau;
make_equilibrium(f,c,t,v,SN2);
out2 = t.mu;





#### Run a season
#include("make_params_s.jl");
#P_C_Sn = eye(P_C_n);
#Fish_xy,Cons_xy,Cons_H = make_season(Fish_xy,Cons_xy,Cons_H,P_C_Sn);

#### Test different social networks at equilibrium
make_params()
initialize()
SN = eye(P_C_n);
make_equilibrium(Fish_xy,Fish_cl,Fclust_xy,Cons_xy,Cons_H,SN);

include("make_params_e.jl");
SN = ones(P_C_n,P_C_n);
make_equilibrium(Fish_xy,Cons_xy,Cons_H,SN);


#### Save
npzwrite("./Data/Data_fish.npy", Fish_xy)
npzwrite("./Data/Data_fisher_xy.npy", Cons_xy)
npzwrite("./Data/Data_fisher_H.npy", Cons_H)
npzwrite("./Data/Data_fclust.npy", Fclust_xy)

#### Useful Julia code
###### Grid fish density
#ii = int(round(Fish_xy[:,:,1]));
#jj,mm,kk = hist2d(ii,Grd_x[1:2:end],Grd_y[1:2:end]);

