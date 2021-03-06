

###### PARAMETERS / CONSTANTS for running the model to equilibrium
function init_equilibrium()
#### Initialize fish
#! all we need are fish location (x,y)
#! the school location (x,y)
#! the index of the school each fish is attached to (ind)
#! the indices of the fish in each school
#! and the KDTREE

@set_constants PRM

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
school_pop=transpose(make_pop(PF_frac,PF_n))
for school=1:PS_n-1
    school_pop=[school_pop, transpose(make_pop(PF_frac,PF_n)) ]
end

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
cons_H   = zeros(Float64,PC_n,length(PF_val)); # catch per species at time t
cons_s   = randn(PC_n);cons_s[cons_s.<0]=-1;cons_s[cons_s.>0]=1;cons_s=int(cons_s)
cons_mi  = int(ones(PC_n));
cons_v   = ones(PC_n) .* PC_v;

#### Initialize social network
cons_sn  = zeros(PC_n,PC_n) #.*max( eps(),PC_lambda); 

#Quantities that have to be measured during simulations
cons_measure = ["Ts"=>zeros(Float64,PC_n), "Tv"=>zeros(Float64,PC_n) ,
    "ts"=>zeros(Float64,PC_n) , "ns"=>zeros(Float64,PC_n), 
    "distance"=>zeros(Float64,PC_n),"turns" =>zeros(Float64,PC_n)]  
cons_series=Dict()

##### Initialize
school = School(school_xy,school_fish,school_pop);
fish = Fish(fish_fx,fish_fs);
cons = Fishers(cons_xy,cons_Ni,cons_target,cons_Nd,cons_dx,
               cons_H,cons_s,cons_mi,cons_sn,cons_v,
               cons_measure,cons_series)
OUT  = Output(fish_fx,cons_xy,school_xy,school_pop,cons_H,cons_mi);



#==== Events & Flags ======#
 #Events list what happened to whom (fish, fisher or school) 
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


#Flags: Switches that control the behaviour of the simulation
 
 FLAGS=(ASCIIString=>Union(Int,Bool,ASCIIString))[
   "benichou"=>true, #  can detect fish only at rest (in search mode)
    "rtree"=>true,  # Use r-tree to find nearest neighbor among fish
    "save"=>false,  # Print out positions of all fishers and fish to make movies
    "spying"=>false,  # unidirectional communciation allowed
    "random_partition"=>false, #every run, partition fishers into random cliques
    "implicit_fish"=>true,  # instead of modelling discrete fish, schools are disks
    "decision"=>false,  # Elaborate decision process (WIP)
    "renewal"=>false,  # Fish are renewed over time
    "measure_frac"=>false,  # measure_frac: measure fraction of time spent in school
    "measure_fij"=>false,  # measure_fij: measure fraction of time spent close to other fishers (SLOW, requires measure_frac)
    "measure_H"=>true,  # measure_H: measure catch rate (including variation over time)
    "measure_flux"=>false,  # measure_flux: measure flux between fisher states
    "stop"=>2,  # Stop condition
    "stopamount"=>convert(Int,PF_n*PS_n.*PC_n*20), # If the stop condition requires a parameter
    "IFQ"=>false,  # Disable fishers once their catch reaches stopamount
    ]


#==== Fishtree =====#
 #One tree per school
if FLAGS["rtree"] && !FLAGS["implicit_fish"]
    fishtree=Fishtree([fnc_makefishtree(i,school,fish) for i=1:PS_n] )
else
    #fake fishtree when deactivated or not needed (no explicit fish)
    fishtree = Fishtree([])
end

return school,fish,cons,fishtree,EVENTS,FLAGS,OUT
end

#============ Initializing Network ============#
function init_network(cons,FLAGS,cliq=[])
    cons_sn=cons.SN
    @set_constants PRM
    #Cliques
    scliq=floor(PC_n/PC_ncliq) #clique size, first approx
    if length(cliq)==0
        #If no clique list is given, divide fishers into
        #equal cliques (number of cliques given by PC_nqcliq)
        cliq=cell(PC_ncliq)
        for c=1:PC_ncliq
            cliq[c]=[(c-1)*scliq+ n  for n=1:scliq ]
        end
        i=1
        n= scliq  * PC_ncliq +1
        while n<PC_n   #distribute remaining fishers
            cliq[i]=[cliq[i], n ]
            n+=1
            i+=1
        end
    end
    for c in cliq
        for i in c
            for j in c
                #Spying flag: Wolf and sheep experiment
                if !FLAGS["spying"] || j< PC_n/2 && i>PC_n/2
                    #Either not spying, or spying and target is in 1st half and spy is in 2nd half
                    cons_sn[i,j] = max( eps(),PC_lambda)
                end
            end
        end
    end
    #Diagonal
    for j = 1:PC_n; cons_sn[j,j] = 1; end;
    cons.SN=cons_sn
end
