	
#### Add modules
using PyPlot
using Types, Constants, NPZ

#### Add functions and routines
include("sub_functions.jl");
include("sub_init.jl");
include("sub_routines.jl");
include("Experiments.jl");

####### Run model for one season
fish,cons,OUT = init_equilibrium();
SN = ones(PC_n,PC_n) .* .01;
for j = 1:PC_n; SN[j,j] = 1; end;

@time make_season(fish,cons,SN,1);

npzwrite("./Data/Data_fish.npy", OUT.fish_xy)
npzwrite("./Data/Data_fishers.npy", OUT.cons_xy)
npzwrite("./Data/Data_clusters.npy", OUT.schl_xy)
npzwrite("./Data/Data_harvest.npy", OUT.cons_H)





###### Simple problem
@time CPUE,Tau = sim_simple();
npzwrite("./Data/Data_simple.npz", ["x"=>1, "CPUE"=>CPUE, "Tau"=>Tau])
x = linspace(0,1,10);
y = squeeze(mean(CPUE,2),2); 
yerr = squeeze(std(CPUE,2),2);
yerr = yerr ./ y[1]; y = y ./ y[1];
errorbar(x,y,xerr=0,yerr=yerr)


###### Optimization problem
#@time (MU, S2) = sim_optimize();
#npzwrite("./Data/Data_opt_cpue_mu.npy", MU)
#npzwrite("./Data/Data_opt_cpue_s2.npy", S2)

###### Evolutionary problem
#@time SN,CPUE_mu,CPUE_s2 = sim_evolve();
#npzwrite("./Data/Data_evo_SN.npy", OUT_SN)
#npzwrite("./Data/Data_evo_cpue_mu.npy", OUT_CPUE_mu)
#npzwrite("./Data/Data_evo_cpue_s2.npy", OUT_CPUE_s2)
##
