#Push module load path to julia pathenv
#push!(LOAD_PATH, "C:/Users/theplankt/Documents/Github/AgentBasedModel/")

using Devectorize #implement devectorizing optimizations
using Types, Constants, Rect

include("sub_functions.jl");
include("sub_init.jl");
include("sub_routines.jl");
include("Experiments.jl");
println("Libraries loaded: working:")

#@profile CPUE,Tau = sim_simple()
#println(CPUE)
#println(Tau)
