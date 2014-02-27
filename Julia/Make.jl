
#### Add modules
using NPZ, Distance, Types, Constants

#### Add functions (i.e. things the model uses)
include("sub_functions.jl");
include("sub_init.jl");

#### Add routines (i.e. the model)
include("sub_routines.jl");

### Initialize
fish,cons,vars,tau,OUT = init_equilibrium();

### make social network
SN = eye(PC_n);
#SN = ones(PC_n,PC_n);

#### Run model
make_equilibrium(fish,cons,tau,vars,SN,1);

#### Save for plotting
npzwrite("./Data/Data_fish.npy", OUT.fish_xy)
npzwrite("./Data/Data_fishers.npy", OUT.cons_xy)
npzwrite("./Data/Data_clusters.npy", OUT.clus_xy)
npzwrite("./Data/Data_harvest.npy", OUT.cons_H)


### Optimization problem
function sim_optimize()
ad = linspace(0,1,10);
TAU_MU = zeros(length(ad));
TAU_S2 = zeros(length(ad));
for i = 1:length(ad)

    ## modulate social network
    SN = ones(PC_n,PC_n) .* ad[i];
    for j = 1:PC_n; SN[j,j] = 1.; end;

    ## Initialize
    fish,cons,vars,tau = init_equilibrium();

    ## run model
    make_equilibrium(fish,cons,tau,vars,SN,0);

    ## calculate average encounter rate
    TAU_MU[i] = mean(tau.mu);
    TAU_S2[i] = mean(tau.s2);

end
return TAU
end

### Run experiments
@time TAU = sim_optimize();




#### Useful Julia code
###### Grid fish density
#ii = int(round(Fish_xy[:,:,1]));
#jj,mm,kk = hist2d(ii,Grd_x[1:2:end],Grd_y[1:2:end]);

#(Tau_n,Tau_t,Tau_s,Tau_mu,Tau_dmu,
#   Dmin,DDx,DDy,ANG,VR,RN,JJ,KK,
#   Fish_xy,Fish_ci,Fish_cl,
#   Cons_xy,Cons_H) = init_equilibrium();
#
#fish = Fish(copy(Fish_xy),copy(Fish_ci),copy(Fish_cl));
#cons = Fishers(copy(Cons_xy),copy(Cons_H));
#vars=Vars(copy(Dmin),copy(DDx),copy(DDy),copy(ANG),copy(VR),copy(RN),copy(JJ),copy(KK));
#tau = Tau(copy(Tau_n),copy(Tau_t),copy(Tau_s),copy(Tau_mu),copy(Tau_dmu));

#npzwrite("./Data/Data_fclust.npy", fish.cl)
#npzwrite("./Data/Data_fisher_H.npy", cons.H)


