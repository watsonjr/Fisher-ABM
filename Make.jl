	
#### Add modules
using PyPlot, NPZ,
@everywhere using Types, Constants, Devectorize

#### Add functions and routines
@everywhere include("sub_functions.jl");
@everywhere include("sub_init.jl");
@everywhere include("sub_routines.jl");
@everywhere include("Experiments.jl");
println("Libraries loaded: working:")

#### Initialize
fish,cons,SN,OUT = init_equilibrium();


####### SIMULATION EXPERIMENTS #######
###### run for one season
@time make_season(fish,cons,SN,1);
npzwrite("./Data/Data_fish.npy", OUT.fish_xy)
npzwrite("./Data/Data_fishers.npy", OUT.cons_xy)
npzwrite("./Data/Data_clusters.npy", OUT.schl_xy)
npzwrite("./Data/Data_harvest.npy", OUT.cons_H)

###### Simple problem (raise SN entirely)
@time CPUE,Tau = sim_simple();
npzwrite("./Data/Data_simple.npz", ["x"=>1, "CPUE"=>CPUE, "Tau"=>Tau])

###### Optimize SN for the fleet catch
@time (CPUE,Social_network) = sim_fleet();

###### Optimize an Individual's social network
@time (CPUE,Social_network) = sim_individual()
npzwrite("./Data/Data_individual.npz", ["x"=>1, "CPUE"=>CPUE, "SN"=>Social_network])
####### END  #######


