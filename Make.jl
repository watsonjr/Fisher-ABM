	
#### Add modules
using PyPlot, NPZ, Devectorize
using Types, Constants

#### Add functions and routines
include("sub_init.jl");
include("sub_functions.jl");
include("sub_routines.jl");
include("Experiments.jl");


##########################!!!!! SIMULATION EXPERIMENTS !!!!!!!!################
###### run for one season
global PC_rp = 0.99; # choose random change in walk
fish,cons,OUT = init_equilibrium();
@time Ts,Tv = make_season(fish,cons,0);
npzwrite("./Data/Data_fish.npy", OUT.fish_xy)
npzwrite("./Data/Data_fishers.npy", OUT.cons_xy)
npzwrite("./Data/Data_clusters.npy", OUT.schl_xy)
npzwrite("./Data/Data_harvest.npy", OUT.cons_H)



##################### 1 FISH 1 FISHER
###### Test performance as a function of random walk
RP = linspace(0.9,1.,20);
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


###### Test performance as a function of C_f and F_sig
SIG = linspace(0,GRD_mx,10);
FF  = linspace(0,GRD_mx,10);
for i = 1:length(Sig)
	for j = 1:length(FF)
		global PF_sig = SIG[i];
		global PC_f   = FF[j];
		fish,cons,OUT = init_equilibrium();
		Ts[i] = make_season(fish,cons,0);
		print(i,"\n")
	end
end

##################### 1 FISH 2 FISHERS
######## Test effect of friendship, on Tau_s and CPUE, as Tau_l is varied


######## Test effect of friendship, on Tau_s and CPUE, as F_n is varied






##################### 1 FISH N FISHERS (do 5, 10, 20 fishers)
######## Stupid Optimal Social Network test (just keep adding links)

######## Evolving social network - find the ESS 








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


