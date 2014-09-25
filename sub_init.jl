###### PRAMETERS / CONSTANTS for running the model to equilibrium
function init_equilibrium()

#### Initialize fish
#! all we need are fish location (x,y)
#! the school location (x,y)
#! the index of the school each fish is attached to (ind)
#! the indices of the fish in each school
#! and the KDTREE

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


#### Initialize school


##### Initialize fishers
#! all we need are the locations of each fisher (x,y)
#! the index of the nearest fish (from KDTree)
#! the index of the target fish (possibly acquired through friend)
#! the distance to that target fish
#! the unit vector pointing at that fish (random at start)
#! the harvest count (1/0)
#! the social strategy (make/break friendships)
#! an index of whether, in the absense of fish, they steam or search
#! the social network
cons_xy  = rand(Float64,PC_n,2).*repmat([GRD_mx GRD_mx],PC_n,1);
cons_Ni  = zeros(Int,PC_n,1);
cons_target  = zeros(Int,PC_n,1);
cons_Nd  = Array(Float64,PC_n,1)
cons_Nd[:]=NaN
DXY      = [randn(PC_n) randn(PC_n)];
cons_dx  = DXY ./ sqrt(DXY[:,1].^2 + DXY[:,2].^2);
cons_H   = zeros(Float64,PC_n); # catch at time t
cons_s   = randn(PC_n);cons_s[cons_s.<0]=-1;cons_s[cons_s.>0]=1;cons_s=int(cons_s)
cons_mi  = int(ones(PC_n));
cons_Ts  = zeros(PC_n);
cons_Tv  = zeros(PC_n);
cons_ts  = zeros(PC_n);
cons_ns  = zeros(PC_n);
cons_sn  = ones(PC_n,PC_n) .*max( eps(),PC_lambda); 
for j = 1:PC_n; cons_sn[j,j] = 1; end;
cons_v   = ones(PC_n) .* PC_v;


##### Initialize
school = School(school_xy,school_fish);
fish = Fish(fish_fx,fish_fs);
cons = Fishers(cons_xy,cons_Ni,cons_target,cons_Nd,cons_dx,
			   cons_H,cons_s,cons_mi,cons_sn,cons_Ts,cons_Tv,
			   cons_ts,cons_ns,cons_v)
OUT  = Output(fish_fx,cons_xy,school_xy,cons_H);



#==== Events & Flags ======#
 #Events list what happened to whom (fish, fisher or school) 
 #within the last timestep to know where updates are needed
 #NB: specific functions control specific events, try not
 #to change them everywhere.
 
 #Events for fish: captured
 #Events for fisher: captor, targeting (fish found but outside harvesting distance)
 #		new neighbor (located another fish), new target (changed targeted fish)
 #Events for school: jumped
 EVENTS=(ASCIIString=>Set{Int})["captured"=>Set{Int}(), "captor"=>Set{Int}(),
 	"new_target"=>Set{Int}(),"new_neighbor"=>Set{Int}(),
 	"targeting"=>Set{Int}(),"jumped"=>Set{Int}()]


#Flags: Switches that control the behaviour of the simulation
 # benichou: can detect fish only at rest
 # rtree: Use r-tree to find nearest neighbor among fish
 FLAGS=(ASCIIString=>Bool)["benichou"=>true,"rtree"=>true,"save"=>false]


#==== Fishtree =====#
 #One tree per school
 fishtree=Fishtree([fnc_makefishtree(i,school,fish) for i=1:PS_n] )

return school,fish,cons,fishtree,EVENTS,FLAGS,OUT
end

#### REFS
# 1 - http://math.stackexchange.com/questions/52869/numerical-approximation-of-levy-flight/52870#52870
# 2 - http://www.johndcook.com/standard_deviation.html
