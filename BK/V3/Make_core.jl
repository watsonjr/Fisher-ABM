#### Add modules
using Types
using NPZ, Devectorize
using PyCall
import Iterators
@pyimport rtree.index as pyrtree #R-tree python module
unshift!(PyVector(pyimport("sys")["path"]), "") #include local folder in search
pypartition=pyimport(:random_partition) #Random partition for experiment

####!! Add functions and routines
include("Utilities.jl");
include("Constants.jl");
include("sub_init.jl");
include("sub_functions.jl");
include("sub_routines.jl");
include("Experiments_core.jl");

####!! Switches for various experiments below
do_1fisher_1school_a = false #expected search time 1 fisher vs C_rp
do_1fisher_1school_b = false #expected search time 1 fisher vs Cf, Fsigma 
do_2fisher_1school = true #VOI over many tau_h, tau_l
do_mfisher_nschool_cliq_a=false  #optimal clique size for 1 tau_h, 1 tau_l
do_mfisher_nschool_cliq_b=false  #optimal clique size over many tau_h, tau_l

####!! Basic timing/profiling test for a single run
if timingtest
    println("First run")
    do_timingtest(true)
    println("Then profiling (now that everything is compiled)")
    @profile do_timingtest(false)
    logfile=open("log.dat","w")
    Profile.print(logfile,format=:flat)
    close(logfile)
end

####!! Experiments for benchmarks paper
##(function definitions are in Experiments.jl)
if do_1fisher_1school_a
    reinit_parameters()
    do_1fisher_1school_a()
end

if do_1fisher_1school_b
    reinit_parameters()
    do_1fisher_1school_b()
end

if do_2fisher_1school
    reinit_parameters()
    do_2fisher_1school()
end

if do_mfisher_nschool_cliq_a
    reinit_parameters()
    do_mfisher_nschool_cliq_a()
end

if do_mfisher_nschool_cliq_b
    reinit_parameters()
    do_mfisher_nschool_cliq_b()
end
