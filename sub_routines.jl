
#### Run a season
function make_season(fish,cons,ST)
 
 #! RUN
 #! while difference in estimated Tau_s is greater that 1%
 #! and the minimum number of schools visited is 10 
 dTs = ones(PC_n);
 while maximum(dTs) > 0.001 || minimum(cons.ns) < 500
 
    ## Distances
    #D,Dx,Dy,cons.MI = fnc_distance(fish.fx,cons.x,cons.MI);
    cons.Ni = fnc_fishfinder(fish.fx,fish.sx,fish.fs,cons.x,GRD_mx2,PC_f);
    #(cons.Ni,cons.Dmin) = fnc_distance_3(fish.fx,cons.x,PC_f);
 
 	## Update steam(MI=0)/search(MI=1) switch
 	cons.MI = fnc_steam(cons.MI)
 
    ## Contact network from probabilistic social network
    CN = fnc_contact(cons.SN)
 
 	## Gather Information
 	#! return nearest distance, updated heading for nearest fish,
 	#! index of nearest fish, harvest success/failure index,
 	(cons.Dmin,cons.DXY,JJ,KK,cons.V) = fnc_information(cons.DXY,cons.Ni,
 	   								fish.fx,cons.x,cons.MI,CN);
 
    ## Harvest
    #! update te cumulative harvest and fish locations
    (cons.H,fish.fx) = fnc_harvest(KK,JJ,cons.H,fish.fx);
 
    ## Move
    #! update positions
    (fish.fx,fish.sx,cons.x) = fnc_move(fish.sx,fish.fx,fish.fs,
     									 cons.x,cons.Dmin,cons.DXY,cons.V);
 
	## Estimate expected time searching for a school
	(cons.Ts,cons.Tv,cons.ts,cons.ns,dTs) = fnc_tau(KK,cons.Ts,
											cons.Tv,cons.ts,cons.ns,dTs);

 	## Save
     if ST == 1
         OUT.fish_xy = cat(3,OUT.fish_xy,fish.fx);
         OUT.cons_xy = cat(3,OUT.cons_xy,cons.x);
         OUT.schl_xy = cat(3,OUT.schl_xy,fish.sx);
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

return cons.Ts, cons.Tv
end


#### Run a season
function make_trip(fish,cons,ST)
 
 #while min of cumulative harvest is less than all fish in the region
 while minimum(cons.H) .< (PF_n*PS_n.*5)
 	#! stops when every fisherman has caught the # fish in the system
 
     ## Distances
     #D,Dx,Dy,cons.MI = fnc_distance(fish.fx,cons.x,cons.MI);
     cons.Ni = fnc_fishfinder(fish.fx,fish.sx,fish.fs,cons.x,GRD_mx2,PC_f);
     #(cons.Ni,cons.Dmin) = fnc_distance_3(fish.fx,cons.x,PC_f);
 
 	 ## Update steam/search switch
 	 cons.MI = fnc_steam(cons.MI)
 
     ## Contact network from probabilistic social network
     CN = fnc_contact(cons.SN)
 
 	 ## Gather Information
 	 #! return nearest distance, updated heading for nearest fish,
 	 #! index of nearest fish, harvest success/failure index,
 	 (cons.Dmin,cons.DXY,JJ,KK) = fnc_information(cons.DXY,cons.Ni,
 										fish.fx,cons.x,cons.MI,CN);
 
     ## Harvest
     #! cons. => CC in function scope
     #! update te cumulative harvest and fish locations
     (cons.H,fish.fx) = fnc_harvest(KK,JJ,cons.H,fish.fx);
 
     ## Move
     #! update positions
     (fish.fx,fish.sx,cons.x) = fnc_move(fish.sx,fish.fx,fish.fs,
     									 cons.x,cons.Dmin,cons.DXY);
 
 	## Storage for plotting
     if ST == 1
         OUT.fish_xy = cat(3,OUT.fish_xy,fish.fx);
         OUT.cons_xy = cat(3,OUT.cons_xy,cons.x);
         OUT.schl_xy = cat(3,OUT.schl_xy,fish.sx);
         OUT.cons_H  = cat(2,OUT.cons_H,cons.H);
     end
 end

end

