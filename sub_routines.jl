### Run multiple seasons 

function make_ensemble(nens,FLAGS,cliq=[])
    #Runs a number "nens" of seasons, keeping the same flags
    measure=Dict()
    series=Dict()
    measure["Hrate_all"]=Float64[]
    measure["Hdist_all"]=Float64[]

    #In case of clique partitions, prepare dicts for catch per clique size
    by_size=Dict{Int,Float64}() #by clique size
    nb_size=Dict{Int,Float64}() #number of cliques of given size

    print("\n")
    for i in 1:nens
        print("|")

        ## If fishers are randomly partitioned into cliques, create partition
        if FLAGS["random_partition"]
            part=pycall(pypartition["random_partition"],PyAny,PRM.PC_n) #random partition of PC_n into integers
            cliq=cell(length(part)) #composition of the cliques
            idx=1
            for p in 1:length(part)
                cliq[p]=idx:idx+part[p]-1
                idx+=part[p]
            end
        end
        
        #If there are cliques, store which fisher is in which clique
        if length(cliq)>0
            belong=zeros(Float64,PRM.PC_n) #clique to which each fisher belongs
            for p in 1:length(cliq)
                for fisher in cliq[p]
                    belong[fisher]=p
                end
            end
        end

        # ============== Make run ==================================

        school,fish,cons,fishtree,EVENTS,ignoredflags,OUT = init_equilibrium();
        init_network(cons,FLAGS,cliq)
        make_season(school,fish,cons,fishtree,EVENTS,FLAGS);

        # ==========================================================

        #compute catch by clique and by clique size
        if length(cliq)>0
            cons.measure["Hclique"]=zeros(Float64,PRM.PC_n) #by clique
            for p in 1:length(cliq)
                by_clique=mean([cons.measure["Hrate"][fisher] for fisher in cliq[p]])
                s=length(cliq[p])
                if in(s,keys(by_size))
                    by_size[s]+=by_clique
                    nb_size[s]+=1
                else
                    by_size[s]=by_clique
                    nb_size[s]=1             
                end
                for f in cliq[p]
                    cons.measure["Hclique"][f]=by_clique
                end
            end
        end
        
        #All the quantities in cons.measure are averaged out
        for k in keys(cons.measure)
            if ! in(k,keys(measure))
                measure[k]=cons.measure[k]        
            else
                measure[k].+=cons.measure[k]
            end
        end
        append!(measure["Hrate_all"],cons.measure["Hrate"][:,1])
        append!(measure["Hdist_all"],cons.measure["Hdist"][:,1])
        
        
        #All the time series in cons.series are concatenated
        for k in keys(cons.series)
            if ! in(k,keys(measure)) || length(series[k])==0
                series[k] = cons.series[k]
            elseif length(cons.series[k])!=0
                series[k] = [series[k],cons.series[k]]
            end
        end
    end
    measure["Hrate_std"]=std(measure["Hrate_all"])*nens
    measure["Hdist_std"]=std(measure["Hdist_all"])*nens
    delete!(measure,"Hrate_all")
    delete!(measure,"Hdist_all")
    for k in keys(measure) #Rescale measured quantities
        measure[k]/=nens
    end
    if length(cliq)>0
        measure["H_by_clique_size"]=zeros(Float64,PRM.PC_n)
        for s in keys(by_size)
            measure["H_by_clique_size"][s]=by_size[s]/nb_size[s]
        end 
        #print("\n")
        #println(measure["H_by_clique_size"])
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


