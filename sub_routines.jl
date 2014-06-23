
#### Run a season
function make_trip(fish,cons,SN,ST)

#while min of cumulative harvest is less than all fish in the region
while minimum(cons.cs) .< (PF_n*PS_n)
            #! stops when every fisherman has caught the # fish in the system

    ## Distances
    #!calculate distances between fish and fishermen? + search/steam switch
    D,Dx,Dy,cons.MI = fnc_distance(fish.fx,cons.x,cons.MI);

    ## Contact network from probabilistic social network
    #!randomly decide if entities contact each other depending on current friendship ?
    fnc_contact(SN,cons.CN); # updates cons.CN

    ## Individual actions
    for i = 1:PC_n

        ## gather information
        #! return nearest distance, updated heading for nearest fish,
        #! index of nearest fish, harvest success/failure index,
        (cons.Dmin[i],cons.DXY[i,:],cons.JJ[i],cons.KK[i]) =
            fnc_information(D,Dx,Dy,cons.DXY[i,:],cons.MI[i],cons.CN,i);

    end

    ## Harvest
    #! cons. => CC in function scope
    #! update te cumulative harvest and fish locations
    (cons.H,cons.cs,fish.fx) = fnc_harvest(cons.KK,cons.JJ,cons.H,cons.cs,fish.fx);

    ## Move
    #! update positions
    (fish.fx,fish.sx,cons.x,cons.Dist) = fnc_move(fish.sx,fish.fx,fish.fs,
    									 cons.x,cons.DXY,cons.Dist);

	## Storage for plotting
    if ST == 1
        OUT.fish_xy = cat(3,OUT.fish_xy,fish.fx);
        OUT.cons_xy = cat(3,OUT.cons_xy,cons.x);
        OUT.schl_xy = cat(3,OUT.schl_xy,fish.sx);
        OUT.cons_H  = cat(2,OUT.cons_H,cons.H);
    end
end
return
end

