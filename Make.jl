



#### Add modules
#using PyPlot
using Types#, Constants
using NPZ, Devectorize
using PyCall
@pyimport rtree.index as pyrtree



#### Add functions and routines
include("Constants.jl");
include("sub_init.jl");
include("sub_functions.jl");
include("sub_routines.jl");
include("Experiments.jl");



#### Switches for various experiments below
timingtest=false
firstpass=false
fig2a=false
fig2b=false
fig3=true



#### Basic timing/profiling test for a single run

if timingtest
    #do_timingtest()
    @profile do_timingtest()
    Profile.print(format=:flat)
end


#### Experiments for benchmarks paper

if firstpass
    do_first_passage()
end

if fig2a
    do_fig2a()
end

if fig2b
    do_fig2b()
end
if fig3
    do_fig3()
end


#################### OLDER EXPERIMENTS #############################

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
