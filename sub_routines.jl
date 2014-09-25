
#### Run a season
function make_season(school,fish,cons,ST,stopflag=2)
 
 #! RUN
 #! while difference in estimated Tau_s is greater that 1%
 #! and the minimum number of schools visited is 10 
 dTs = ones(PC_n);
 
 #=== Possible loop conditions ==============#
 
 cond1=(fish,cons,dTs)-> (minimum(cons.H) .< (PF_n*PS_n.*5))
 #! while min of cumulative harvest is less than all fish in the region
 #! stops when every fisherman has caught the # fish in the system
 
 cond2=(fish,cons,dTs)->(maximum(dTs) > .001 || minimum(cons.ns) < 500)
 #! while difference in estimated Tau_s is greater that 0.1%
 #! and the minimum number of schools visited is <500 
 
 whilecond=[cond1,cond2][stopflag]

 
 
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
 FLAGS=(ASCIIString=>Bool)["benichou"=>true,"rtree"=>true]

 #==== Fishtree =====#
 #One tree per school
 fishtree=Fishtree([fnc_makefishtree(i,school,fish) for i=1:PS_n] )
 
 args=[school,fish,cons,fishtree,EVENTS,FLAGS] # arguments of most functions
 turns=0
 while whilecond(fish,cons,dTs)
 	turns+=1
    ## Distances
    #D,Dx,Dy,cons.MI = fnc_distance(fish.fx,cons.x,cons.MI);
    fnc_fishfinder(GRD_mx2,PC_f,args...);
    #(cons.Ni,cons.Dmin) = fnc_distance_3(fish.fx,cons.x,PC_f);
 
 	## Update steam(MI=0)/search(MI=1) switch
 	fnc_steam(args...);
 
    ## Contact network from probabilistic social network
    ## This network includes only people who have new information
    ## and are willing to share it (random trial)
    CN = fnc_contact(args...)
 
 	## Gather Information
 	#! return nearest distance, updated heading for nearest fish,
 	#! index of nearest fish, harvest success/failure index,
 	fnc_information(CN,args...);
 
    ## Harvest
    #! update te cumulative harvest and fish locations
    fnc_harvest(args...);
 
    ## Move
    #! update positions
    fnc_move(args...);
 
	if stopflag==2
		## Estimate expected time searching for a school
		fnc_tau(dTs,cons,EVENTS);
	end
	
 	## Save
     if ST == 1
         OUT.fish_xy = cat(3,OUT.fish_xy,fish.fx);
         OUT.cons_xy = cat(3,OUT.cons_xy,cons.x);
         OUT.schl_xy = cat(3,OUT.schl_xy,school.x);
         OUT.cons_H  = cat(2,OUT.cons_H,cons.H);
     end
 end
# Save
if ST == 1
	npzwrite("./Data/Data_fish.npy", OUT.fish_xy)
	npzwrite("./Data/Data_fishers.npy", OUT.cons_xy)
	npzwrite("./Data/Data_clusters.npy", OUT.schl_xy)
	npzwrite("./Data/Data_harvest.npy", OUT.cons_H)
end

#return cons.Ts, cons.Tv # expectation and variance in time between schools
end


