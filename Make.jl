	
#### Add modules
using PyPlot, NPZ, Devectorize
using Types, Constants

#### Add functions and routines
include("sub_init.jl");
include("sub_functions.jl");
include("sub_routines.jl");
include("Experiments.jl");


####### SIMULATION EXPERIMENTS #######
###### run for one season
#@time make_trip(fish,cons,1);
global PC_rp = 1.;
fish,cons,OUT = init_equilibrium();
@time Ts = make_season(fish,cons,1);
npzwrite("./Data/Data_fish.npy", OUT.fish_xy)
npzwrite("./Data/Data_fishers.npy", OUT.cons_xy)
npzwrite("./Data/Data_clusters.npy", OUT.schl_xy)
npzwrite("./Data/Data_harvest.npy", OUT.cons_H)


###### Test different behaviors
RP = linspace(0,1,10);
Ts = cell(size(RP));
for i = 1:length(RP)
	global PC_rp = RP[i];
	fish,cons,OUT = init_equilibrium();
	Ts[i] = make_season(fish,cons,0);
	print(i,"\n")
end



TS = zeros(size(RP))
for i = 1:length(RP)
	TS[i] = Ts[i][1]
end
npzwrite("./Data/Data_Ts.npy", TS)



####### Simple problem (raise SN entirely)
#@time CPUE,Tau = sim_simple();
#npzwrite("./Data/Data_simple.npz", ["x"=>1, "CPUE"=>CPUE, "Tau"=>Tau])
#
####### Optimize SN for the fleet catch
#@time (CPUE,Social_network) = sim_fleet();
#
####### Optimize an Individual's social network
#@time (CPUE,Social_network) = sim_individual()
#npzwrite("./Data/Data_individual.npz", ["x"=>1, "CPUE"=>CPUE, "SN"=>Social_network])
######## END  #######


