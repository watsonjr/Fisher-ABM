#### Add modules
using Types
using NPZ, Devectorize
using PyCall
import Iterators
@pyimport rtree.index as pyrtree #R-tree python module
unshift!(PyVector(pyimport("sys")["path"]), "") #include local folder in search
pypartition=pyimport(:random_partition) #Random partition for experiment

####!! Add functions and routines
include("Constants.jl");
include("sub_functions.jl");
include("sub_init.jl");
include("sub_routines.jl");
include("Experiments_core.jl");

####!! Switches for various experiments below
one_fisher_one_school_a = false #expected search time 1 fisher vs C_rp
one_fisher_one_school_b = false #expected search time 1 fisher vs Cf, Fsigma 
two_fisher_one_school = true #VOI over many tau_h, tau_l
mfisher_nschool_cliq_a=false  #optimal clique size for 1 tau_h, 1 tau_l
mfisher_nschool_cliq_b=false  #optimal clique size over many tau_h, tau_l

####!! Basic timing/profiling test for a single run
println("First run")
do_timingtest(true)
println("Then profiling (now that everything is compiled)")
@profile do_timingtest(true)
logfile=open("log.dat","w")
Profile.print(logfile,format=:flat)
close(logfile)

####!! Experiments for benchmarks paper
##(function definitions are in Experiments.jl)
if one_fisher_one_school_a
    reinit_parameters()
    do_1fisher_1school_a()
end

if one_fisher_one_school_b
    reinit_parameters()
    do_1fisher_1school_b()
end

if two_fisher_one_school
    reinit_parameters()
    tic(); do_2fisher_1school(); toc()
end

if mfisher_nschool_cliq_a
    reinit_parameters()
    do_mfisher_nschool_cliq_a()
end

if mfisher_nschool_cliq_b
    reinit_parameters()
    do_mfisher_nschool_cliq_b()
end
