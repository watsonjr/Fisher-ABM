
#### Run a season
function make_trip(fish,cons,ST)
 
 #while min of cumulative harvest is less than all fish in the region
 while minimum(cons.H) .< (PF_n*PS_n)
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

