

########======== Make An Ensemble of Seasons
function make_ensemble(school,fish,cons,fishtree,EVENTS,FLAGS,OUT=None)
    nens=FLAGS["ensemble"]
    if nens==1
         make_season(school,fish,cons,fishtree,EVENTS,FLAGS,OUT)
    else
        all=school,fish,cons,fishtree,EVENTS,FLAGS,OUT
        collected=Array(Any,nens)
        for ens=1:nens
            result=deepcopy(all)
            make_season(result...)
            collected[ens]=result[1:3]
        end
        #Take the average over ensembles for measurements in cons.measure
        for i in cons.measure
            cons.measure[i[1]]=mean( [collected[ens][3].measure[i[1]] for ens=1:nens] )
        end
        #Take the average for all fields in school, 
		#fish or cons that are arrays of floats
        for i in 1:3
            @simd for f in names(all[i])
                if issubtype( typeof(all[i].(f)),Array{Float64})
                    test=mean( [getfield(collected[ens][i],f) for ens=1:nens] )
                    @inbounds all[i].(f)=test
                end
            end
        end
    end
end



########============= Make One Season
function make_season(school,fish,cons,fishtree,EVENTS,FLAGS,OUT=None)
 #======= Preamble =======#
 @set_constants PRM
 dHs = ones(PC_n); dTs = ones(PC_n); turns=0;
 args=[school,fish,cons,fishtree,EVENTS,FLAGS] # arguments of most functions
 
 #======= Loop Conditions =======#
 if FLAGS["whilecond"] == 1
 	# stop if every fisher has encountered tonnes of schools
 	whilecond = (fish,cons,dTs,dHs,EVENTS)-> (minimum(cons.H)<(PF_n.*PS_n.*30))
 elseif FLAGS["whilecond"] == 2
 	# stop if expected harvest rate has converged
 	whilecond = (fish,cons,dTs,dHs,EVENTS)->(maximum(dHs) > .01 ||
        	minimum(cons.measure["ns"]) < 500)
 end

 #======= RUN MODEL TO EQUILIBRIUM ========#
 while whilecond(fish,cons,dTs,dHs,EVENTS)
 	## Ticker
    turns+=1

    ## Distances
    #D,Dx,Dy,cons.MI = fnc_distance(fish.fx,cons.x,cons.MI);
    fnc_fishfinder(args...);
    #(cons.Ni,cons.Dmin) = fnc_distance_3(fish.fx,cons.x,PC_f);
 
    ## Update steam(MI=0)/search(MI=1) switch
    fnc_steam(args...);
 
    ## Contact network from probabilistic social network
    ## This network includes only people who have new information
    ## and are willing to share it (random trial)
    CN = fnc_contact(args...);
 
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
 
    ## Estimate expected time searching for a school
    fnc_tau(dTs,dHs,cons,EVENTS,FLAGS,turns);

    ## Save for a movie
    if FLAGS["save"] == true
         OUT.fish_xy  = cat(3,OUT.fish_xy,fish.fx);
         OUT.cons_xy  = cat(3,OUT.cons_xy,cons.x);
         OUT.schl_xy  = cat(3,OUT.schl_xy,school.x);
         OUT.schl_pop = cat(2,OUT.schl_pop,school.pop);
         OUT.cons_H   = cat(2,OUT.cons_H,cons.H);
    end
     
end

# Save for making movie
if FLAGS["save"] == true
    npzwrite("./Data/Data_fish.npy", OUT.fish_xy)
    npzwrite("./Data/Data_fishers.npy", OUT.cons_xy)
    npzwrite("./Data/Data_clusters.npy", OUT.schl_xy)
    npzwrite("./Data/Data_cluspop.npy", OUT.schl_pop)
    npzwrite("./Data/Data_harvest.npy", OUT.cons_H)
    fout=open("./Data/Data_figs.dat", "w")
    for f in names(PRM)
        write(fout,"$f    $(getfield(PRM,f))\n" )
    end
    close(fout)
end

return turns
end


