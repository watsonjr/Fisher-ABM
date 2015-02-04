

###### PARAMETERS / CONSTANTS for running the model to equilibrium
function init_equilibrium()
#### Initialize fish
#! all we need are fish location (x,y)
#! the school location (x,y)
#! the index of the school each fish is attached to (ind)
#! the indices of the fish in each school
#! and the KDTREE
@set_constants PRM

#### Initialize random location of fish and fishers
fish_fx = zeros(Float64,PF_n*PS_n,2); # location of fish
school_xy=rand(Float64,PS_n,2).*repmat([GRD_mx GRD_mx],PS_n,1);
a = zeros(PS_n,PF_n);
for i = 1:PS_n; a[i,:] = ones(PF_n) .* i; end
fish_fs = a'[:]; # fish school index
if PS_n==1
    school_fish=cat(2,1:PF_n,1:PF_n ) #create false second school so that school_fish[:,1] refers to first school
else
    school_fish=cat(2,[find(fish_fs.==school) for school=1:PS_n]...); #fishes in each school
end
for i = 1:PF_n*PS_n; fish_fx[i,:,1] = mod(school_xy[fish_fs[i],:,1] +
            randn(1,2).*PF_sig, [GRD_mx GRD_mx]); end
school_pop=[PF_n for  school=1:PS_n] # init schools


##### Initialize fishers
#! all we need are the locations of each fisher (x,y)
#! the index of the nearest fish (from KDTree)
#! the index of the target fish (possibly acquired through friend)
#! the distance to that target fish
#! the unit vector pointing at that fish (random at start)
#! the harvest count (1/0)
#! the social network
cons_xy  = rand(Float64,PC_n,2).*repmat([GRD_mx GRD_mx],PC_n,1);
cons_Ni  = zeros(Int,PC_n,1); # index of nearest fish
cons_target = zeros(Int,PC_n,1); # index of target fish
cons_Dmin = Array(Float64,PC_n,1); cons_Dmin[:] = NaN # distance to fish
DXY      = [randn(PC_n) randn(PC_n)]; # initial heading
cons_dx  = DXY ./ sqrt(DXY[:,1].^2 + DXY[:,2].^2);
cons_H   = zeros(Float64,PC_n); # catch at time t
cons_mi  = int(ones(PC_n)); # steam(1)/search(0) switch
cons_sn  = zeros(PC_n,PC_n) # social network
cons_v   = ones(PC_n) .* PC_v; # speed

##### Quantities that have to be measured during simulations
#! Ts: expected time between schools
#! Tv: variance in time between schools
#! ts: needed calculate the running mean
#! ns: number of schools found
#! distance: total distance traveled
#! Hrate: number of fish per unit time
#! Hdist: number of fish per distance
cons_measure = ["Ts"=>zeros(Float64,PC_n), "Tv"=>zeros(Float64,PC_n) ,
    "ts"=>zeros(Float64,PC_n) , "ns"=>zeros(Float64,PC_n), 
    "distance"=>zeros(Float64,PC_n), "Hrate"=>zeros(Float64,PC_n),
    "Hdist"=>zeros(Float64,PC_n) ]  


##### Initialize
school = School(school_xy,school_fish,school_pop);
fish   = Fish(fish_fx,fish_fs);
cons   = Fishers(cons_xy,cons_Ni,cons_target,cons_Dmin,cons_dx,
               cons_H,cons_mi,cons_sn,cons_v,cons_measure)
OUT    = Output(fish_fx,cons_xy,school_xy,school_pop,cons_H);

##### EVENTS:list what happened to whom (fish, fisher or school) 
#within the last timestep to know where updates are needed
#NB: specific functions control specific events, try not
#to change them everywhere.

#Events for fish: captured
#Events for fisher: captor, targeting (fish found but outside harvesting distance)
#        new neighbor (located another fish), new target (changed targeted fish)
#          bound (targeting a fish given by someone else)
#Events for school: jumped
EVENTS=(ASCIIString=>Set{Int})["captured"=>Set{Int}(), "captor"=>Set{Int}(),
     "new_target"=>Set{Int}(),"new_neighbor"=>Set{Int}(),
     "targeting"=>Set{Int}(),"jumped"=>Set{Int}(),
     "in_contact"=>Set{Int}(), "spying"=>Set{Int}(),
     "left_school"=>Set{Int}(),"found_school"=>Set{Int}(),"bound"=>Set{Int}() ]


##### FLAGS: Switches that control the behaviour of the simulation
# benichou: can detect fish only at rest (in search mode)
# rtree: Use r-tree to find nearest neighbor among fish
# save: Print out positions of all fishers and fish to make movies
# implicit_fish: instead of modelling discrete fish, schools are disks
# spying: unidirectional communciation allowed
# measure_frac: measure fraction of time spent in school
# measure_fij: measure fraction of time spent close to other fishers (SLOW, requires measure_frac)
# measure_H: measure catch rate
FLAGS=(ASCIIString=>Any)["benichou"=>true,"rtree"=>true,"save"=>false,
    "spying"=>false,"implicit_fish"=>true,
    "measure_frac"=>false,"measure_fij"=>false,"measure_H"=>true,
    "ensemble"=>1,"whilecond"=>1]

##### Fishtree
#One tree per school
if FLAGS["rtree"] && !FLAGS["implicit_fish"]
    fishtree=Fishtree([fnc_makefishtree(i,school,fish) for i=1:PS_n] )
else
    #fake fishtree when deactivated or not needed (no explicit fish)
    fishtree = Fishtree([])
end

return school,fish,cons,fishtree,EVENTS,FLAGS,OUT
end
