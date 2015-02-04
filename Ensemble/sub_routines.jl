
#### Run a season
function make_ensemble(school,fish,cons,fishtree,EVENTS,FLAGS,stopflag=2,OUT=None)
    nens=FLAGS["ensemble"]
    if nens==1
         make_season(school,fish,cons,fishtree,EVENTS,FLAGS,stopflag,OUT)
    else
        all=school,fish,cons,fishtree,EVENTS,FLAGS,stopflag,OUT
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
        #Take the average for all fields in school, fish or cons that are arrays of floats
        for i in 1:3
            for f in names(all[i])
                if issubtype( typeof(all[i].(f)),Array{Float64})
                    test=mean( [getfield(collected[ens][i],f) for ens=1:nens] )
                    all[i].(f)=test
                end
            end
        end
    end
end

function make_season(school,fish,cons,fishtree,EVENTS,FLAGS,stopflag=2,OUT=None)
 
 @set_constants PRM
 #println("PC_rp")
 #! RUN
 #===== Possible loop conditions ======#
 cond1=(fish,cons,dTs,dHs,EVENTS)-> (minimum(cons.H) < (PF_n*PS_n.*.5))
 #! while min of cumulative harvest is less than all fish in the region
 #! stops when every fisherman has caught the # fish in the
 
 dTs = ones(PC_n);
 cond2=(fish,cons,dTs,dHs,EVENTS)->(maximum(dTs) > .001 || 
 	minimum(cons.measure["ns"]) < 500)
 #! while difference in estimated Tau_s is greater that 0.1%
 #! and the minimum number of schools visited is <500 
 
 dHs = ones(PC_n);
 cond3=(fish,cons,dTs,dHs,EVENTS)->(maximum(dHs) > .01 || 
 	minimum(cons.measure["ns"]) < 500)
 #! Same as cond2 with estimated catchrate
 
 #( maximum(cons.H) < 1)
 cond4=(fish,cons,dTs,dHs,EVENTS)->isempty(EVENTS["found_school"])
 #! stop as soon as a school has been found#old=catch has been 
 #! made (useful for first-passage time)
 
 #! Collect conditions and arguments
 whilecond=[cond1,cond2,cond3,cond4][stopflag]
 args=[school,fish,cons,fishtree,EVENTS,FLAGS] # arguments of most functions
 turns=0
 
 
 #======= RUN MODEL TO EQUILIBRIUM ========#
 while whilecond(fish,cons,dTs,dHs,EVENTS)
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
     if FLAGS["save"] == 1
         OUT.fish_xy  = cat(3,OUT.fish_xy,fish.fx);
         OUT.cons_xy  = cat(3,OUT.cons_xy,cons.x);
         OUT.schl_xy  = cat(3,OUT.schl_xy,school.x);
         OUT.schl_pop = cat(2,OUT.schl_pop,school.pop);
         OUT.cons_H   = cat(2,OUT.cons_H,cons.H);
     end
     
end

# Save for making movie
if FLAGS["save"] == 1
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

#if FLAGS["measure_H"]
#    cons.measure["Hdist"]=cons.H ./ cons.measure["distance"]
#end

#return cons.Ts, cons.Tv # expectation and variance in time between schools
return turns
end


