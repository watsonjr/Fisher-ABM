
######## Functions for ABM



######## Basic Utility Functions ############################

#### Spatial difference function (Vectorized, currently unused)
#! x1 = first x,y location
#! x2 = second x,y location
#! dx,dy = difference in x,y accounting for periodic boundary
function fnc_difference_vec(x1,x2)
    @set_constants PRM
    # difference
    dx = x1[:,1] .- x2[:,1]
    dy = x1[:,2] .- x2[:,2]
    # periodic boundary
    j = (abs(dy).>GRD_mx2) + (abs(dx).>GRD_mx2);
    dx[j.>1] = -sign(dx[j.>1]).*(GRD_mx .- abs(dx[j.>1]));
    dy[j.>1] = -sign(dy[j.>1]).*(GRD_mx .- abs(dy[j.>1]));

    return dx,dy
end

#### Spatial difference and distance functions (Devectorized)
function fnc_difference(pos1,pos2)
    GRD_mx2= PRM.GRD_mx2
    dpos=pos1-pos2
    df=div(dpos,GRD_mx2)
    #println("test",GRD_mx2," ",df, " ",((dpos)%GRD_mx2-df*GRD_mx2))
    return ((dpos)%GRD_mx2-df*GRD_mx2)
end

function fnc_dist(pos1,pos2)
    GRD_mx2= PRM.GRD_mx2
    dpos=pos1-pos2
    return hypot( (dpos)%GRD_mx2-div(dpos,GRD_mx2)*GRD_mx2... )
end

#Disable a fisher entirely (stops searching around, harvesting, communicating)
function fnc_disable_fisher(i,cons,EVENTS)
    cons.MI[i]=-1 #disabled mode
    cons.V[i]=0
    cons.target[i]=0
    cons.Ni[i]=0
    for k in keys(EVENTS) #Remove disabled fisher from all events
        delete!(EVENTS[k],i)
    end
end

#### OPTIMAL TURN RATE
#! Analytical expression of the optimal turn rate probability
#! as a function of other global parameters
function fnc_optimal_PCrp()
    @set_constants PRM
    v=PC_v
    a=PC_f+PF_sig#min(PF_sig,PF_n*PC_f/(2*pi) )
    b=(a+GRD_mx/sqrt(PS_n) )/2.
    if log(b/a)>.5
        tau2=max(1.01,(a/v * sqrt( log(b/a) -1./2.)))
    else
        println("$a $b $(PS_n) $(PF_sig) $(PS_n*pi*PF_sig^2/GRD_mx^2) ")
        return 0.001
    end
    return 1.-1./tau2
end

#### Benichou estimate for tau_s of a single fisher
function fnc_taus1()
    @set_constants PRM
    v=PC_v
    a=PC_f+PF_sig
    b=(a+GRD_mx/sqrt(PS_n) )/2.
    t2=1./(1.-PC_rp) #avg time of straight flight
    t1=1./(PC_rp) #avg time of search
    xy=sqrt(2)/v/t2
    x=a*xy
    y=b*xy
    result1=xy^1/a*t1*(b^2-a^2)^2*besseli(0,x)/besseli(1,x)
    result2=xy^2 *t1*(4*b^4*log(b/a) + (b^2-a^2)*(a^2-3*b^2+8/xy^2))
    result3= result1 + result2/4.
    return (t1+t2)/(2*t1*b^2) * result3
end

### tau_h for a single fisher and tau_l (harvesting time and landscape change time)
function taucalc()
    @set_constants PRM
    tauh=PF_n/PC_q
    taul=1./PS_p
    return tauh,taul
end

#### Get identity of school associated with a target
function get_school(tgt,fish,FLAGS)
    if FLAGS["implicit_fish"]
        return tgt
    else
        return fish.fs[tgt]
    end
end

######## ABM MEASUREMENT FUNCTIONS ############################
# Functions that are used in simulations to measure some quantities 


# Initialize all arrays that are used to contain measured quantities
# They are activated by flags
function init_measurements(cons,FLAGS,EVENTS)
    @set_constants PRM
 cons.measure["math"]=zeros(4) #f s b bbar
 if FLAGS["measure_frac"]
    cons.measure["f1"]=zeros(PC_n)
    cons.measure["f2"]=zeros(PC_n)
    cons.measure["fij"]=zeros(PC_n)
    EVENTS["finders"]=Set{Int}()
 end
 if FLAGS["measure_H"]
    cons.series["Hflux"]=[0]
 end
 
 if FLAGS["measure_flux"]
    for i in ("found_school","left_school","lost_school", "bound")
        cons.series[i]=[]
    end
 end
 cons.measure["dTs"] = ones(PC_n); #Difference in tau_s
 cons.measure["dHs"] = ones(PC_n); #Difference in catch rate
end

# At the end of a run, wrap up the measurements
function wrap_measurements(cons,FLAGS,OUT)
    @set_constants PRM
    # Export data for movie
    if FLAGS["save"] == 1
        npzwrite("./Data/Data_fish.npy", OUT.fish_xy)
        npzwrite("./Data/Data_fishers.npy", OUT.cons_xy)
        npzwrite("./Data/Data_clusters.npy", OUT.schl_xy)
        npzwrite("./Data/Data_cluspop.npy", OUT.schl_pop)
        npzwrite("./Data/Data_harvest.npy", OUT.cons_H)
        npzwrite("./Data/Data_MI.npy", OUT.cons_MI)
        fout=open("./Data/Data_figs.dat", "w")
        for f in names(PRM)
            write(fout,"$f    $(getfield(PRM,f))\n" )
        end
        close(fout)
    end

    turns=cons.measure["turn"]

    if FLAGS["measure_frac"]
        cons.measure["f1"]/=turns
        cons.measure["f2"]/=turns
        cons.measure["fij"]/=turns
        cons.measure["math"]./=turns
    end

    if FLAGS["measure_H"]
        Htot=sum(cons.H.*transpose(PF_val) ,2)
        cons.measure["Htot"]=Htot
        cons.measure["Hdist"]=Htot./ cons.measure["distance"]
        cons.measure["Hrate"]=Htot./cons.measure["turns"]
        cons.series["Hflux"]=cons.series["Hflux"][2:end]-cons.series["Hflux"][1:end-1]
    end
end

#### HAUL TIME
#! running time between hauls
#! and estimate the running mean time between schools
#! and estimate the difference in this running mean 
#! which is the switch for the while loop
function fnc_tau(cons,EVENTS,FLAGS)
    @set_constants PRM
    Ts,Tv,ts,ns,dTs,dHs=cons.measure["Ts"],cons.measure["Tv"],cons.measure["ts"],cons.measure["ns"],cons.measure["dTs"],cons.measure["dHs"]
    turns=cons.measure["turn"]

    for I=EVENTS["left_school"]
        ts[I]=0
        if FLAGS["measure_H"]
            cons.measure["dHs"][I]=-sum(cons.H.*transpose(PF_val))/turns +10. #store catch rate when leaving school (1)
        end
    end
    
    for I=EVENTS["found_school"]
        ns[I] += 1; # update school counter
        Ts_old = Ts[I]; # current mean
        Tv_old = Tv[I]; # current variance

        Ts[I] = Ts_old + ((ts[I]-Ts_old)/ns[I]) # run mean
        Tv[I] = Tv_old + ((ts[I]-Ts_old)*(ts[I]-Ts[I])) # run variance
        
        dTs[I] = abs(Ts[I]-Ts_old)/Ts[I]; # fractional change in mean
        ts[I] = 1; # reset how long it took to find school
        if FLAGS["measure_H"]
            if sum(cons.H[I,:])>1
                dHs[I]=abs(1+ (dHs[I]-10.)/sum(cons.H[I,:].*PF_val )*turns) 
                #compare (1) to catch rate when reaching next school
            end
        end      
    end
    if FLAGS["measure_H"]
        cons.series["Hflux"]=[cons.series["Hflux"],sum(cons.H.*transpose(PF_val)) ]
    end
    for I = 1:PC_n
        ts[I]+=1
    end

    return
end



######## ABM CORE FUNCTIONS ############################
# Functions that are used in every simulation (part of the routine)



#### MAKE FISHTREE 
#! (active only if FLAGS["rtree"]==true)
#! Everytime a school jumps, or at initalization
#! create the corresponding kd-tree
#! for fast nearest-neighbor lookup
function fnc_makefishtree(i,school,fish)
    @set_constants PRM
    tree=pyrtree.Index()
    for f = school.fish[:,i] #for all fish in the school
        fx=fish.fx[f,:]
        tree[:insert](f,hcat(fx,fx+PC_h )') 
        #insert an element to the tree that has the index f of the fish as label
        #and the position of the fish as a location (rtrees require quadrangles so
        #the position is doubled)
    end
    return tree
end


##### When an event happens (found school, etc.), append agent to appropriate list
## And if events are recorded, note the current state (number of searching, bound, etc. fishers)
function record_event(event,i,EVENTS,FLAGS,cons)
    if event in keys(EVENTS)
        push!(EVENTS[event],i)
    end
    if FLAGS["measure_flux"] && "states" in keys(cons.series)
        tauh,taul=taucalc()
        tmp=[cons.series["states"][end,:] tauh taul ]  #Record current states
        if length(cons.series[event])>0
            cons.series[event]=[cons.series[event],tmp] 
        else
            cons.series[event]=tmp
        end
    end
end

#### Typical distance between fishermen
function measure_dist(cons)
    @set_constants PRM
    D=[]
    if PC_n==1
        return 0
    end
    for i=2:PC_n
        for j=1:i
            if length(D)>0
                D=[D fnc_dist( cons.x[i,:],cons.x[j,:])]
            else
                D=fnc_dist( cons.x[i,:],cons.x[j,:])
            end
        end
    end
    return [mean(D) std(D)]
end

#### FISH FINDER / DISTANCE
#!calculate distances between fishermen and nearest target if any within scope
#!(targets can be fish in the base model, or schools in the "implicit_fish" model)
function fnc_fishfinder(school,fish,cons,fishtree,EVENTS,FLAGS)
    @set_constants PRM
    empty!(EVENTS["new_neighbor"])
    empty!(EVENTS["left_school"])
    empty!(EVENTS["found_school"])
    
    #fout=open("./Data/EVTS.dat", "a")

    grd=GRD_mx/2

    Fx,Sx,Si,Cx=fish.fx,school.x,fish.fs,cons.x

    # First, find nearby schools
    II = cell(PC_n); # index of schools close to each fisher
    for i = 1:PC_n
        II[i] =[]
        if FLAGS["benichou"] &&   cons.MI[i]!=1 
            ## Benichou mode: look for neighbors only if fisher currently in search mode
            continue
        end
        
        if FLAGS["implicit_fish"]
            #Neighbors are schools, find nearest neighbor
            newk=0
            Dmin=NaN
            for j=1:PS_n
                D=fnc_dist( Cx[i,:],Sx[j,:])
                if D<(PC_f+PF_sig)  &&  (isnan(Dmin) || D< Dmin)
                    newk=j
                    Dmin=D
                end
            end
            if newk!=cons.Ni[i]
                if newk==0
                    record_event("left_school",i,EVENTS,FLAGS,cons)

                else
                    push!(EVENTS["new_neighbor"],i)
                    record_event("found_school",i,EVENTS,FLAGS,cons)
                end
                cons.Ni[i]=newk
            end
        else
            #Neighbors are fish
            #look ford schools that are potentially close enough
            #that some fish may be within reach
            for j=1:PS_n
                D=fnc_dist(Cx[i,:],Sx[j,:])
                if D.<((2.5.*PF_sig).+PC_f)
                   II[i]=[II[i],j]
                end
            end
        end
    end
    #close(fout)

    if FLAGS["implicit_fish"]
        #Neighbors are schools, already obtained above, no need to go further
        return
    end


    # Output - index of nearest fish for all fishers (if any)
    Ni = cons.Ni;
    ts=cons.measure["ts"];
    # find for each fisher
    for i = 1:PC_n

        if isempty(II[i]) == 0 # if there is a school nearby fisher i

            # if so, get index k of fish in all schools, near fisher i
            newk=0
            Dmin=NaN #minimal distance to another fish
            for j = 1:length(II[i])
                #Find nearest neighbor in school j
                if FLAGS["rtree"]
                    k=fishtree.trees[j][:nearest]( hcat(Cx[i,:],Cx[i,:] )' )[:next]()
                    D=fnc_dist(Cx[i,:],Fx[k,:]) 
                    if D<PC_f
                        Dmin=D
                        newk=k
                    end
                    treek=newk
                else
                    for k = school.fish[:,II[i][j]]
                        D=fnc_dist(Cx[i,:],Fx[k,:]) 
                        if D<PC_f && (isnan(Dmin) || D<Dmin)
                            newk=k
                            Dmin=D
                        end
                    end
                end
            end
            if newk!=Ni[i]
                if newk!=0
                    if Ni[i]==0 && ts[i]>4*min(PC_f,PF_sig)/PC_v
                        push!(EVENTS["found_school"],i)
                    end
                    push!(EVENTS["new_neighbor"],i)
                else
                    push!(EVENTS["left_school"],i)
                end
                Ni[i]=newk
            end
        else 
            if Ni[i]!=0
                Ni[i] = 0;
                push!(EVENTS["left_school"],i)
            end
        end
    end

    #cons.Ni=Ni # Index of nearest fish for each fisher (in fishfinder)
    #return Ni
end




#### Search/steam switch
#! if a fisher can't see any fish
#! they can either spin around or move in a straightish line
#! they switch between this probabilistically
function fnc_steam(school,fish,cons,fishtree,EVENTS,FLAGS)
    @set_constants PRM
    MI=cons.MI
    V=cons.V
    for i =  1:PC_n
        if cons.target[i]!=0 || MI[i]==-1 #Fish found or fisher disabled
            continue
        end
        
        if MI[i] == 1 # if searching
            if true#rand() < PC_rp # maybe switch to steaming
                MI[i] = 0
                V[i]=PC_v
            end
        elseif MI[i] == 0 # if steaming
            if rand() < PC_rp #(1-PC_rp) # maybe switch to searching
                MI[i] = 1
                V[i]=0#PC_v/3.
                #Watch out for the diffusive velocity in search mode
                #(arbitrarily set to PC_v/3 here because that fits
                #the benichou curves with k=infinity best, but changing
                #this doesn't change the optimal value of PC_rp anyway)
            end
        end
    end
    return MI
end


#### CONTACT NETWORK from social network
#! Iterates through social network adjacency matrix, generates 2 random numbers.
#! If random numbers are less than adjacency measure (friendship) they both
#! contact; else both do not contact. Return the network 2D array.
function fnc_contact(school,fish,cons,fishtree,EVENTS,FLAGS)
    @set_constants PRM
    SN=cons.SN
    CN = zeros(Int,PC_n,PC_n)
    Cx=cons.x
    
    empty!(EVENTS["in_contact"])
    probing= EVENTS["new_neighbor"]  #ask to any friend who has a new neighbor
    
    for j in probing
        CN[j,j]=1
        for i = 1:PC_n
            if i==j || cons.MI[i]==-1
                continue
            end
            f1 = SN[i,j];
            f2 = SN[j,i];
            RN = rand(2,1);
           
            if FLAGS["spying"] 
                dx = fnc_dist(Cx[i,:],Cx[j,:])  #distance between fishermen
                #Directional exchanges are possible within spying radius PC_spy
                if RN[1] < f1 && dx < PC_spy
                    CN[i,j] = 1;
                    push!(EVENTS["spying"],i)
                    push!(EVENTS["in_contact"],i)
                end
            else
                #Sharing only
                if RN[1] < f1 && RN[2] < f2
                    CN[i,j] = 1;
                    CN[j,i] = 1;
                    push!(EVENTS["in_contact"],i)
                    push!(EVENTS["in_contact"],j)
                end
            end
        end
    end
    return CN
end
 
#### DECISION PROCESS for the fishermen (UNUSED FOR NOW)
#! As a function of the signal sent by a source (self or
#! other fisher) and the distance to the target,
#! determine the expected remaining number of fish by
#! the time one reaches the target, then compare it
#! to current expected catch (either from another target
#! or from random motion, knowing average proba of finding school)
#! Two parameters: memory, allowing to integrate a signal from the
#! same source over time; knowledge of terrain
function fnc_decision(id,Fx,CN,school,fish,cons,fishtree,EVENTS,FLAGS)
    current=cons.Dmin[id]
    dxy,Ni,Cx,MI=cons.DXY,cons.Ni,cons.x,cons.MI
    Dmin=NaN
    if in(id,EVENTS["in_contact"]) #If you have just talked with someone
        # get the vector of people with whom you are currently in contact
        J  = find(CN[id,:].==1); # index of friends
        Jn = length(J); # number of friends 
        # calculate distances to all targets you have info on
        DD = fill(NaN,Jn) # distance to your target and friend's
        Dx = fill(NaN,Jn) # dx
        Dy = fill(NaN,Jn) # dy
        TGT=fill(NaN,Jn) #target
        for i = 1:Jn # for each friend
            ii = J[i]; # get his/her index
            jj = Ni[ii]; #get their fish
            TGT[i]=jj
            if jj != 0 && !isnan(Fx[jj,1]) # if they see a fish
                (dx,dy) = fnc_difference(Fx[jj,:],Cx[id,:]);
                DD[i] = hypot(dx,dy); # and distance
                Dx[i] = dx; Dy[i] = dy;
            end
        end
    else
        #Only your own fish
        if  Ni[id]!=0 && isnan(current) 
            J=[id,]
            Jn=1
            jj=Ni[id]
            (dx,dy) = fnc_difference(Fx[jj,:],Cx[id,:]);
            Dx=[dx,]; Dy=[dy,];
            DD=[hypot(dx,dy),]
            TGT=[id,]
        else
            Jn=0
            return cons.target[id],cons.DXY[id],cons.Dmin[id]
        end
    end
    Dmin = NaN
    Dxy=  cons.DXY[id,:]

    target=0
    if !FLAGS["decision"]
        #Basic decision process: go to closest
        if Jn>0
            ii= indmin(DD) #nearest
            jj = Ni[J[ii]]
            Dxy = [Dx[ii] Dy[ii]] ./ norm([Dx[ii] Dy[ii]]);  # calculate unit vector DXY to nearest fish
            Dmin = DD[ii]; # shortest distance to a target
        end
        if ( !isnan(Dmin)  && (isnan(current) || Dmin<current ) )
            target=jj
        elseif  in(id,EVENTS["targeting"]) && !isnan(current) 
            #Stay with current target
            target=cons.target[id]
            Dxy=cons.DXY[id,:]
            Dmin=cons.Dmin[id]
        end
    else
        #Advanced decision process
        taus=fnc_taus1()
    end
    return target,Dxy,Dmin
    
end

#### INFORMATION and fisher direction
#! Iterate through the people you share information with and get the locations
#! of the fish within their view. Calculate the unit vector to the nearest fish.
#! Otherwise roam around randomly according to a self-correlated walk that 
#! approximates search behavior according to a Levy walk
#! Return the minimum distance, updated heading, JJ index of nearest fish
function fnc_information(CN,school,fish,cons,fishtree,EVENTS,FLAGS)
    @set_constants PRM
    #In this function we manage targeting events
    
    #fout=open("./Data/EVTS.dat", "a")
    
    if FLAGS["implicit_fish"]
        #Neighbors are schools rather than fish
        Fx=school.x
    else
        Fx=fish.fx
    end

    dxy,Ni,Cx,MI=cons.DXY,cons.Ni,cons.x,cons.MI
    DMIN,DXY=cons.Dmin,cons.DXY # shortest distance & unit vector
    V     = cons.V # speed

    #Fishers for whom previous target has been canceled:
    #fish captured or part of school that jumped
    for id=EVENTS["targeting"]
        tgt=cons.target[id]
        if in(tgt, EVENTS["captured"]) || in(get_school(tgt,fish,FLAGS), EVENTS["jumped"]) 
            cons.Dmin[id]=NaN
            if tgt==cons.Ni[id]
                cons.Ni[id]=0
            end
            cons.target[id]=0
            #if FLAGS["measure_frac"]
            #    delete!(EVENTS["finders"],id)
            #end
            delete!(EVENTS["targeting"],id)
            
            if FLAGS["measure_flux"]
                if in(id,EVENTS["bound"])
                    record_event("lost_school",id,EVENTS,FLAGS,cons)
                end
            end
            delete!(EVENTS["bound"],id)
            V[id]=PC_v #Restore speed to normal value
        end
    end
    
    for id = 1:PC_n
        Dxy=DXY[id,:]
        Dmin=DMIN[id]
        target,Dxy,Dmin=fnc_decision(id,Fx,CN,school,fish,cons,fishtree,EVENTS,FLAGS)
        # Decide which target to choose
        if target!=cons.target[id]
            # if I see anything better than current target (if any)
            if in(id,EVENTS["targeting"])
                #In case the fisher was previously bound
                delete!(EVENTS["bound"],id)
            end
            
            if target!=cons.Ni[id]
                #If the target is not own neighbour
                record_event("bound",id,EVENTS,FLAGS,cons)
            elseif FLAGS["measure_frac"]
                push!(EVENTS["finders"],id)
            end
            push!(EVENTS["targeting"],id)
            
            cons.target[id]=target
            cons.V[id]=PC_v
        elseif target==0
            # else I roam around randomly 
            cons.target[id]=0
            delete!(EVENTS["targeting"],id)
            if FLAGS["measure_frac"]
                delete!(EVENTS["finders"],id)
            end

            #Steam or search pattern
            if MI[id] == 0 # # steam
                #angle = atand(dxy[1]/dxy[2]);
                #angle = angle .+ (PC_r*2) * rand() .- PC_r;
                #dx  = cosd(angle); dy  = sind(angle);
                #Dxy = [dx dy];
                #Dxy = dxy[id,:] + (randn(1,2).*PC_r) # bendy walk
                #Dxy = Dxy ./ norm(Dxy)
                Dxy = dxy[id,:]; # straight line
                #V[id] = PC_v # fast speed
            else # search/tumble
                #angle = 360*rand();
                #dx  = cosd(angle); dy  = sind(angle);
                #Dxy = [dx dy];
                #Dxy = dxy[id,:] + (randn(1,2).*PC_r2)
                Dxy = randn(1,2); # random walk
                Dxy = Dxy ./ norm(Dxy)
                #V[id] = PC_v *0 # slow speed
            end
        end
        # slow speed
        #if ! isnan(Dmin) && Dmin<PC_h
        #    V[id] = PC_v / 5
        #end
        # Store
        DXY[id,:] = Dxy
        DMIN[id] = Dmin
    end
    #close(fout)

end


#### HARVEST
function fnc_harvest(school,fish,cons,fishtree,EVENTS,FLAGS);
    @set_constants PRM
    empty!(EVENTS["captor"]) #In this function we manage capture events
    empty!(EVENTS["captured"])
    
    CH,FX,TGT,DMIN=cons.H,fish.fx,cons.target,cons.Dmin
    
    nbcaptors=zeros(PC_n)
    Dharvest=PC_h #max harvest distance
    if FLAGS["implicit_fish"]
        Dhvarest=PF_sig+PC_h
    end
    
    Dmin=cons.Dmin
    
    schools=zeros(Int,PC_n)
    if FLAGS["measure_frac"]
        #Measure time spent in school
        f1=cons.measure["f1"] #Time spent by 1 fisher in a school
        f2=cons.measure["f2"] #Time spent by 2+ fishers in the same school
    end
    if FLAGS["measure_fij"]
        fij=cons.measure["fij"] #Time spent by 2 fishers within distance x
    end

    #For all fishers who have a target
    for i = EVENTS["targeting"]
        if Dmin[i]<Dharvest #If the target is within harvesting distance
            tgt=TGT[i]
            
            if FLAGS["measure_frac"]
                #Add to time spent in school
                schools[i]=get_school(tgt,fish,FLAGS)
            end
            
            if FLAGS["implicit_fish"]
                #If the school is a simple disk with a population variable
                #cons.V[i]=0 #stop ship (maybe not necessary)
                species=school.pop[tgt,:]/sum(school.pop[tgt,:]) #array of current species
                delta=min(school.pop[tgt,:],PC_q*species)  #use probability of catch as depletion rate
                school.pop[tgt,:]-=delta
                cons.H[i,:]+=delta
            else
                #For explicit fish
                #TODO: Deal with mutiple species!
                fx=FX[tgt,:] #position of the target
                if isnan(fx[1]) #if the same fish was captured by someone else
                    continue
                end
                cons.V[i]=0 #The boat stops moving until it catches the fish
                # probabilistic catch
                r =rand()
                if r < PC_q
                    push!(EVENTS["captured"],tgt)
                    school.pop[fish.fs[tgt]]-=1 #School depleted by one
                    if FLAGS["rtree"]
                        #Fish is removed from the tree
                        fishtree.trees[fish.fs[tgt] ][:delete](tgt,hcat(fx,fx+PC_h )') 
                    end
                    FX[tgt,:]=NaN; #Fish disappears
                    push!(EVENTS["captor"],i)
                    nbcaptors[i]+=1
                    cons.V[i] = PC_v
                end
            end
        end
    end
    
    #For explicit fish
    if !FLAGS["implicit_fish"]
        for i=EVENTS["captor"]
            CH[i] += 1./nbcaptors[i]; #Fisher gets proportion of fish
        end
    end

    #If the schools are renewed instead of jumping
    if FLAGS["implicit_fish"] & FLAGS["renewal"]
        school.pop=min(PF_n,school.pop+PS_p.*PF_frac)
    end

    #In case of individual fishing quotas
    #disable fisher if they have reached their quota
    if FLAGS["IFQ"] 
        for i in EVENTS["targeting"]
           if sum(cons.H[i,:].*PF_val)>=FLAGS["stopamount"]
                fnc_disable_fisher(i,cons,EVENTS)
            end
        end
    end

    
    if FLAGS["measure_frac"] ||  FLAGS["measure_flux"]
        #Fraction of time spent by:
        #f1: one fisher in a school
        #f2: two fishers in the same school
        #fij: two fishers within "school size" of each other
        
        #Number of fishers per school
        schfreq=hist(schools,0.5:PS_n+0.5 )[2]
        tmp = [length(EVENTS["finders"]) length(EVENTS["bound"]) length(find(schfreq)) sum(schools.>0)]
        cons.measure
        cons.measure["math"][3]+=tmp[3] #Number of occupied schools
        cons.measure["math"][4]+=tmp[4] #number of people in school
        
        #Number of fishers in a school with k fishers
        fk=hist(schfreq,0.5:PC_n+0.5 )[2]
        
        cons.measure["math"][1]+=tmp[1]
        cons.measure["math"][2]+=tmp[2]
        if !( "states" in keys(cons.series))
            cons.series["states"]=tmp
        else
            cons.series["states"]=[cons.series["states"],tmp]
        end
        D=[cons.measure["turn"] measure_dist(cons)]
        if "dist" in keys(cons.series)
            cons.series["dist"]=[cons.series["dist"],D]
        else
            cons.series["dist"]=D
        end

    end
        
    if FLAGS["measure_frac"] 
        for i=1:PC_n
            if schools[i]>0
                others=schfreq[schools[i]]-1 #Number of other fishers in same school
                if others==0
                    f1[i]+=1
                elseif others>0
                    f2[i]+=1
                end
            end
            if FLAGS["measure_fij"]
                for j=1:PC_n
                    if i!=j && fnc_dist(cons.x[i,:],cons.x[j,:])<PC_f+PF_sig
                        fij[i]+=1
                    end
                end
            end
        end
    end
end



#### MOVE
#! CL <- school.x (school location); FX <- fish location;
#! distance
#! Randomly move the fish schools. 
#! Return the locations of the fish, school locations,
#! updated fisher locations
function fnc_move(school,fish,cons,fishtree,EVENTS,FLAGS)
    @set_constants PRM
    #CL=school.x;FX=fish.fx;FS=fish.fs;CC=cons.x;Dm=cons.Dmin;DXY=cons.DXY;
    CL,FX,CC,Dm,DXY,V=school.x,fish.fx,cons.x,cons.Dmin,cons.DXY,cons.V
    # schools and fish move
    empty!(EVENTS["jumped"])
    for i = 1:PS_n
        if sum(school.pop[i,:])<1 || rand() < PS_p # if no fish left in school or randomly, jump
            if FLAGS["renewal"]
                #Renewal rather than jump
                push!(EVENTS["jumped"],i);
                continue
            end
            j = school.fish[:,i];
            CL_x = rand() * GRD_mx; #center location
            CL_y = rand() * GRD_mx;
            CL[i,:] = [CL_x CL_y];
            school.pop[i,:]=transpose(make_pop(PF_frac,PF_n)) #restore population counter to base level
            push!(EVENTS["jumped"],i);
            if FLAGS["implicit_fish"]
                continue
            end
            F_x  = mod(CL_x.+(randn(PF_n,1)*PF_sig),GRD_mx);
            F_y  = mod(CL_y.+(randn(PF_n,1)*PF_sig),GRD_mx);
            FX[j,:] = [F_x F_y];
            if FLAGS["rtree"] 
                fishtree.trees[i]=fnc_makefishtree(i,school,fish)
            end
        end
    end


    for i=find(x->x!=-1,cons.MI) # Excepting all disabled fishers (MI==-1)
        cons.measure["turns"][i]+=1 #increase time at sea
    end
    
    # fishers move
            # slow down as you approach fish
            #v = Dm; v[v.<PC_h] = PC_h; v[v.>PC_f] = PC_f;
            #v = (v .- PC_h) ./ (PC_f-PC_h);
            #range2 = PC_v - PC_vmn;
            #V = (v*range2) .+ PC_vmn;

     
     for f=EVENTS["targeting"]
         #reduce distance to target
         tgt=cons.target[f]
         dx,dy=fnc_difference(school.x[tgt,:], cons.x[f,:])
         dxy=[dx dy] ./ norm([dx dy])
         if cons.Dmin[f]>cons.V[f]
             cons.Dmin[f]-=cons.V[f]
         else
             #Slow down approach to avoid overshooting the target
             cons.V[f]= cons.Dmin[f]
             cons.Dmin[f]=0
         end
     end

     for f=1:PC_n
         CC_x = mod(CC[f,1] + (DXY[f,1]*V[f]), GRD_mx);
         CC_y = mod(CC[f,2] + (DXY[f,2]*V[f]), GRD_mx);
         cons.x[f,:] =[CC_x CC_y]
         cons.measure["distance"][f]+=V[f]
     end

end


