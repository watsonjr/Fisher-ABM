
#### Add modules
using Types, Constants, Distributions
using NPZ, PyPlot

#### Add functions and routines
include("sub_functions.jl");
include("sub_init.jl");
include("sub_routines.jl");
include("Experiments.jl");

####### Run model for one season
fish,cons,OUT = init_equilibrium();
SN = ones(PC_n,PC_n) .* 0.1;
for j = 1:PC_n; SN[j,j] = 1; end;
make_equilibrium(fish,cons,SN,1);
npzwrite("./Data/Data_fish.npy", OUT.fish_xy)
npzwrite("./Data/Data_fishers.npy", OUT.cons_xy)
npzwrite("./Data/Data_clusters.npy", OUT.clus_xy)
npzwrite("./Data/Data_harvest.npy", OUT.cons_H)


###### Optimization problem
#@time (MU, S2) = sim_optimize();
#npzwrite("./Data/Data_opt_cpue_mu.npy", MU)
#npzwrite("./Data/Data_opt_cpue_s2.npy", S2)
#
#
###### Evolutionary problem
#@time SN,CPUE_mu,CPUE_s2 = sim_evolve();
#npzwrite("./Data/Data_evo_SN.npy", OUT_SN)
#npzwrite("./Data/Data_evo_cpue_mu.npy", OUT_CPUE_mu)
#npzwrite("./Data/Data_evo_cpue_s2.npy", OUT_CPUE_s2)
##
