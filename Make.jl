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

timingtest=false #timing & profiling
firstpass=false #Benichou test: time of first passage
fig2a=false #Fig2a: search time 1 fisher vs C_rp
fig2b=true #Fig2b: search time 1 fisher vs Cf, Fsigma 
fig2c=true #Fig2b: search time 1 fisher vs S_p, C_q 
fig3=true #Fig3: VOI against tau_h/tau_s, tau_l/tau_s
fig4opt=false  #optimal lambda
fig4opt_cliq=false #optimal nb cliques
fig4opt_comp=false #cliques vs lambda
rndcliq=false  #random partition of fishers into cliques
worst=false #like fig3 + iteration over lambda to find min(VOI)

spying=false #spying radius

#### Basic timing/profiling test for a single run

if timingtest
    println("First run")
    do_timingtest(true)
    println("Then profiling (now that everything is compiled)")
    @profile do_timingtest(false)
    logfile=open("log.dat","w")
    Profile.print(logfile,format=:flat)
    close(logfile)
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

if fig2c
    reinit_parameters()
    do_fig2c()
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

if fig4opt_comp
    reinit_parameters()
    do_fig4opt_comp()
end

if rndcliq
    reinit_parameters()
    do_rndcliq()
end

if worst
    reinit_parameters()
    do_worst()
end

if spying
    reinit_parameters()
    do_spying()
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
