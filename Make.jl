#### Add modules
#using PyPlot
using Types
using NPZ, Devectorize
using PyCall
import Iterators
@pyimport rtree.index as pyrtree #R-tree python module
unshift!(PyVector(pyimport("sys")["path"]), "") #include local folder in search
pypartition=pyimport(:random_partition) #Random partition for experiment


#### Add functions and routines
include("Utilities.jl");
include("Constants.jl");
include("sub_init.jl");
include("sub_functions.jl");
include("sub_routines.jl");
include("Experiments.jl");

#### Switches for various experiments below
timingtest=false
firstpass=false
fig2a=true
fig2b=true
fig3=true
fig4opt=true
fig4opt_cliq=true


#### Basic timing/profiling test for a single run

if timingtest
    #do_timingtest()
    @profile do_timingtest()
    Profile.print(format=:flat)
end


#### Experiments for benchmarks paper
##(function definitions are in Experiments.jl)

if firstpass
    reinit_parameters()
    do_first_passage()
    
    reinit_parameters()
    do_fstpass_nschool()
end

if fig2a
    reinit_parameters()
    do_fig2a()
end

if fig2b
    reinit_parameters()
    do_fig2b()
end

if fig3
    reinit_parameters()
    do_fig3()
end

if fig4opt
    reinit_parameters()
    do_fig4opt()
end


if fig4opt_cliq
    reinit_parameters()
    do_fig4opt_cliq()
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
