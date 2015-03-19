### Run multiple seasons 

function make_ensemble(nens,FLAGS)
    #Runs a number "nens" of seasons, keeping the same flags

    #First iteration
    #(Less painful than having to initialize the dictionary "measure" by hand due to types
    school,fish,cons,fishtree,EVENTS,initflags,OUT = init_equilibrium();
    init_network(cons,FLAGS)
    make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
    measure=cons.measure
    series=cons.series
    
    #Next iterations
    for i in 1:nens-1
        school,fish,cons,fishtree,EVENTS,wrongflags,OUT = init_equilibrium();
        init_network(cons,FLAGS)
        make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
        
        #All the quantities in cons.measure are averaged out
        for k in keys(cons.measure)
            measure[k].+=cons.measure[k]
        end
        
        #All the time series in cons.series are concatenated
        for k in keys(cons.series)
            if length(series[k])==0
                series[k] = cons.series[k]
            elseif length(cons.series[k])!=0
                series[k] = [series[k],cons.series[k]]
            end
        end
    end
    for k in keys(measure) #Rescale measured quantities
        measure[k]/=nens
    end
    
    return measure,series
end


#### Run a season
function make_season(school,fish,cons,fishtree,EVENTS,FLAGS,OUT=None)
 
 @set_constants PRM
 #println("PC_rp")
 #! RUN
 #=== Possible loop conditions ==============#
 
 
 cond1=(fish,cons,EVENTS,FLAGS)-> (sum(cons.H[:,1]) < FLAGS["stopamount"] )
 #! while min of cumulative harvest is less than the stop amount
 
 cond2=(fish,cons,EVENTS,FLAGS)->(maximum(cons.measure["dTs"]) > .001 || minimum(cons.measure["ns"]) < 500)
 #! while difference in estimated Tau_s is greater that 0.1%
 #! and the minimum number of schools visited is <500 
 
 cond3=(fish,cons,EVENTS,FLAGS)->(maximum(cons.measure["dHs"]) > .01 || minimum(cons.measure["ns"]) < 500)
 #! Same as cond2 with estimated catchrate
 
 cond4=(fish,cons,EVENTS,FLAGS)->isempty(EVENTS["found_school"]) #( maximum(cons.H) < 1)
 #! stop as soon as a school has been found#old=catch has been made (useful for first-passage time)
 
 cond5=(fish,cons,EVENTS,FLAGS)->(maximum(cons.MI)>-1)
 #! Stop as soon as every fisher has been disabled
 
 whilecond=[cond1,cond2,cond3,cond4,cond5][FLAGS["stop"]]
 
 args=[school,fish,cons,fishtree,EVENTS,FLAGS] # arguments of most functions
 cons.measure["turn"]=0
 
 save_parameters() #Function that writes the parameter values for this run into a file
  
 #Specific measurements for analytical purposes
 init_measurements(cons,FLAGS,EVENTS)
 
 while whilecond(fish,cons,EVENTS,FLAGS)
    cons.measure["turn"]+=1

    ## Distances
    #D,Dx,Dy,cons.MI = fnc_distance(fish.fx,cons.x,cons.MI);
    fnc_fishfinder(args...);
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
 

    ## Estimate expected time searching for a school
    fnc_tau(cons,EVENTS,FLAGS);



    ## Save data for movie
     if FLAGS["save"] == 1
         OUT.fish_xy = cat(3,OUT.fish_xy,fish.fx);
         OUT.cons_xy = cat(3,OUT.cons_xy,cons.x);
         OUT.schl_xy = cat(3,OUT.schl_xy,school.x);
         OUT.schl_pop  = cat(3,OUT.schl_pop,school.pop);
         OUT.cons_H  = cat(3,OUT.cons_H,cons.H);
         OUT.cons_MI = cat(2,OUT.cons_MI,cons.MI);
     end
     
 end

 
 wrap_measurements(cons,FLAGS,OUT)


 return cons.measure["turn"]
end


