#### Add modules
using Types
using NPZ, Devectorize
import Iterators
using PyCall
@pyimport rtree.index as pyrtree #R-tree python module
unshift!(PyVector(pyimport("sys")["path"]), "") #include local folder in search
pypartition=pyimport(:random_partition) #Random partition for experiment

####!! Add functions and routines
include("Constants.jl");
include("sub_functions.jl");
include("sub_init.jl");
include("sub_routines.jl");
include("Experiments.jl");

####!! Switches for various experiments below
timingtest = false #make a movie
one_fisher = false #expected search time 1 fisher vs C_rp
two_fisher = false #VOI over many tau_h, tau_l
nfisher_quantile = false #optimal clique size when VOI high
nfisher_wsls = true #optimal clique size when VOI high

####!! Basic timing/profiling test for a single run
#println("First run")
#do_timingtest(true)
#println("Then profiling (now that everything is compiled)")
#@profile do_timingtest(true)
#logfile=open("log.dat","w")
#Profile.print(logfile,format=:flat)
#close(logfile)

####!! Experiments for benchmarks paper
##(function definitions are in Experiments.jl)
if timingtest
	reinit_parameters()
	do_timingtest(true)
end

if one_fisher
    reinit_parameters()
    do_1fisher()
end

if two_fisher
    reinit_parameters()
    do_2fisher();
end

if nfisher_quantile
    reinit_parameters()
    do_nfisher_quantile(0.9)
end

if nfisher_wsls
    reinit_parameters()
    @time do_nfisher_wsls()
end
