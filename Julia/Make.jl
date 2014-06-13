	
#### Add modules
using PyPlot
using Types, Constants, NPZ

#### Add functions and routines
include("sub_functions.jl");
include("sub_init.jl");
include("sub_routines.jl");
include("Experiments.jl");

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
####### END  #######


###### Quick Plot
x = linspace(0,1,10);
y = squeeze(mean(CPUE,2),2); 
yerr = squeeze(std(CPUE,2),2);
yerr = yerr ./ y[1]; y = y ./ y[1];
errorbar(x,y,xerr=0,yerr=yerr)


