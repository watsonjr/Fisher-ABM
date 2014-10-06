
######## Functions for ABM



######## Basic Utility Functions ############################

#### Spatial difference function (Vectorized, currently unused)
#! x1 = first x,y location
#! x2 = second x,y location
#! dx,dy = difference in x,y accounting for periodic boundary
function fnc_difference_vec(x1,x2)
    # difference
    dx = x1[:,1] .- x2[:,1]
    dy = x1[:,2] .- x2[:,2]
    # periodic boundary
    j = (abs(dy).>GRD_mx2) + (abs(dx).>GRD_mx2);
    dx[j.>1] = -sign(dx[j.>1]).*(GRD_mx .- abs(dx[j.>1]));
    dy[j.>1] = -sign(dy[j.>1]).*(GRD_mx .- abs(dy[j.>1]));

    return dx,dy
end


#### Spatial differnce and distance functions (Devectorized)
function fnc_difference(pos1,pos2)
    dpos=pos1-pos2
    df=div(dpos,GRD_mx2)
    return ((dpos)%GRD_mx2-df*GRD_mx2).*([1 1]-2*df)
end

function fnc_dist(pos1,pos2)
    dpos=pos1-pos2
    return hypot( (dpos)%GRD_mx2-div(dpos,GRD_mx2)*GRD_mx2... )
end


#### OPTIMAL TURN RATE
#! Analytical expression of the optimal turn rate probability
#! as a function of other global parameters
function fnc_optimal_PCrp()
    v=PC_v
    a=PC_f+min(PF_sig,PF_n*PC_f/(2*pi) )
    b=GRD_mx2
    tau2=(a/v * sqrt( log(b/a) -1./2.))
    return 1.-1./tau2
end



######## ABM MEASUREMENT FUNCTIONS ############################
# Functions that are used in some simulations to measure some quantities 

#### HAUL TIME
#! running time between hauls
#! and estimate the running mean time between schools
#! and estimate the difference in this running mean 
#! which is the switch for the while loop
function fnc_tau(dTs,cons,EVENTS)
    Ts,Tv,ts,ns=cons.measure["Ts"],cons.measure["Tv"],cons.measure["ts"],cons.measure["ns"]
    #H=cons.H;Ts=cons.Ts;ts=cons.ts;sn=cons.sn;

    for I=EVENTS["left_school"]
        ts[I]=0
    end
    for I=EVENTS["found_school"]
        #println(ns,ts,Ts)
        ns[I] += 1; # update school counter
        Ts_old = Ts[I]; # current mean
        Tv_old = Tv[I]; # current variance

        Ts[I] = Ts_old + ((ts[I]-Ts_old)/ns[I]) # run mean
        Tv[I] = Tv_old + ((ts[I]-Ts_old)*(ts[I]-Ts[I])) # run variance
        
        dTs[I] = abs(Ts[I]-Ts_old)/Ts[I]; # fractional change in mean
        ts[I] = 1; # reset how long it took to find school
    end
    for I = 1:PC_n
        ts[I]+=1
    end
    return
end



######## ABM CORE FUNCTIONS ############################
# Functions that are used in every simulation (main behavior)



#### MAKE FISHTREE 
#! (active only if FLAGS["rtree"]==true)
#! Everytime a school jumps, or at initalization
#! create the corresponding kd-tree
#! for fast nearest-neighbor lookup
function fnc_makefishtree(i,school,fish)
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




#### FISH FINDER / DISTANCE
#!calculate distances between fishermen and nearest target if any within scope
#!(targets can be fish in the base model, or schools in the "implicit_fish" model)
function fnc_fishfinder(grd,PC_f,school,fish,cons,fishtree,EVENTS,FLAGS)
    empty!(EVENTS["new_neighbor"])
    empty!(EVENTS["left_school"])
    empty!(EVENTS["found_school"])

    Fx,Sx,Si,Cx=fish.fx,school.x,fish.fs,cons.x
    # First, find nearby schols
    II = cell(PC_n); # index of schools close to each fisher
    for i = 1:PC_n
        Dmin=cons.Dmin[i]
        II[i] =[]
        if FLAGS["benichou"] &&   cons.MI[i]!=1 
            ## Benichou mode: look for neighbors only if fisher currently in search mode
            continue
        end
        
        if FLAGS["implicit_fish"]
            #Neighbors are schools, find nearest neighbor
            newk=0
            for j=1:PS_n
                D=fnc_dist( Cx[i,:],Sx[j,:])
                if D<(PC_f+PF_sig)  &&  (isnan(Dmin) || D< Dmin)
                    newk=j
                    Dmin=D
                end
            end
            if newk!=cons.Ni[i]
                if newk==0
                    push!(EVENTS["left_school"],i)
                else
                    push!(EVENTS["new_neighbor"],i)
                    push!(EVENTS["found_school"],i)
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
#! that is, is a fisher can't see any fish
#! they can either spin around or move in a straightish line
#! they switch between this probabilistically
function fnc_steam(school,fish,cons,fishtree,EVENTS,FLAGS)
    MI=cons.MI
    for i =  1:PC_n
        if cons.target[i]!=0
            continue
        end
        #For fishers that don't have a target
        if MI[i] == 1 # if steaming
            if rand() < PC_rp # maybe switch to tumbling
                MI[i] = 0
            end
        elseif MI[i] == 0 # if tumbling
            if rand() < (1-PC_rp) # maybe switch to steaming
                MI[i] = 1
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
    SN=cons.SN
    CN = zeros(Int,PC_n,PC_n)
    Cx=cons.x
    
    empty!(EVENTS["in_contact"])
    probing= EVENTS["new_neighbor"]  #ask to any friend who has a new neighbor
    
    for i = 1:PC_n
        for j in probing
            if i==j
                push!(EVENTS["in_contact"],i)
                CN[i,i]=1
                continue
            end
            f1 = SN[i,j];
            f2 = SN[j,i];
            RN = rand(2,1);
            dx = fnc_dist(Cx[i,:],Cx[j,:])  #distance between fishermen

            if FLAGS["spying"] && dx < PC_spy
                #Directional exchanges are possible within spying radius PC_spy
                if RN[1] < f1
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
function fnc_decision(source,target,school,fish,cons,fishtree,EVENTS,FLAGS)
    return 1
end

#### INFORMATION and fisher direction
#! Iterate through the people you share information with and get the locations
#! of the fish within their view. Calculate the unit vector to the nearest fish.
#! Otherwise roam around randomly according to a self-correlated walk that 
#! approximates search behavior according to a Levy walk
#! Return the minimum distance, updated heading, JJ index of nearest fish
function fnc_information(CN,school,fish,cons,fishtree,EVENTS,FLAGS)
    #In this function we manage targeting events
    
    dxy,Ni,Fx,Cx,MI=cons.DXY,cons.Ni,fish.fx,cons.x,cons.MI
    
    if FLAGS["implicit_fish"]
        #Neighbors are schools rather than fish
        Fx=school.x
    end

    DMIN,DXY=cons.Dmin,cons.DXY # shortest distance & unit vector
    V     = cons.V # speed

    #Fishers for whom previous target has been canceled:
    #fish captured or part of school that jumped
    for id=EVENTS["targeting"]
        tgt=cons.target[id]
        if FLAGS["implicit_fish"]
            school=tgt
        else
            school=fish.fs[tgt]
        end
        if in(tgt, EVENTS["captured"]) || in(tgt, EVENTS["jumped"]) 
            cons.Dmin[id]=NaN
            cons.target[id]=0
            delete!(EVENTS["targeting"],id)
            V[id]=PC_v #Restore speed to normal value
        end
    end
    
    for id = 1:PC_n
        current=cons.Dmin[id]
        Dmin=NaN
        if in(id,EVENTS["in_contact"]) #If you have just talked with someone
            # get the vector of people with whom you are currently in contact
            J  = find(CN[id,:].==1); # index of friends
            Jn = length(J); # number of friends 

            # calculate distances to all targets you have info on
            DD = fill(NaN,PC_n) # distance to your target and friend's
            Dx = fill(NaN,PC_n) # dx
            Dy = fill(NaN,PC_n) # dy
            for i = 1:Jn # for each friend
                ii = J[i]; # get his/her index
                jj = Ni[ii]; #get their fish
                if jj != 0 && !isnan(Fx[jj,1]) # if they see a fish
                    (dx,dy) = fnc_difference(Fx[jj,:],Cx[id,:]);
                    DD[i] = hypot(dx,dy); # and distance
                    Dx[i] = dx; Dy[i] = dy;
                end
            end
            ii= indmin(DD) #nearest
            jj = Ni[J[ii]]
            Dxy = [Dx[ii] Dy[ii]] ./ norm([Dx[ii] Dy[ii]]);  # calculate unit vector DXY to nearest fish
            Dmin = DD[ii]; # shortest distance to a target
        elseif false
            #Only your own fish
            if  Ni[id]!=0 && isnan(current) 
                jj=Ni[id]
                (dx,dy) = fnc_difference(Fx[jj,:],Cx[id,:]);
                Dmin = hypot(dx,dy); # and distance
                Dxy = [dx dy] ./ norm([dx dy]);
            else
                Dmin = NaN
                Dxy=  cons.DXY[id,:]
            end
        end

        # Decide which target to choose
        if  !isnan(Dmin)  && (isnan(current) || Dmin<current ) 
            # if I see anything better than current target (if any)
            push!(EVENTS["targeting"],id)
            cons.target[id]=jj
                
        elseif  in(id,EVENTS["targeting"]) && !isnan(current) 
            #Stay with current target
            Dxy=cons.DXY[id,:]
            Dmin=cons.Dmin[id]
        else 
            # else I roam around randomly 
            cons.target[id]=0
            delete!(EVENTS["targeting"],id)

            #Steam or search pattern
            if MI[id] == 0 # # steam
                #angle = atand(dxy[1]/dxy[2]);
                #angle = angle .+ (PC_r*2) * rand() .- PC_r;
                #dx  = cosd(angle); dy  = sind(angle);
                #Dxy = [dx dy];
                #Dxy = dxy[id,:] + (randn(1,2).*PC_r) # bendy walk
                #Dxy = Dxy ./ norm(Dxy)
                Dxy = dxy[id,:]; # straight line
                V[id] = PC_v # fast speed
            else # tumble
                #angle = 360*rand();
                #dx  = cosd(angle); dy  = sind(angle);
                #Dxy = [dx dy];
                #Dxy = dxy[id,:] + (randn(1,2).*PC_r2)
                Dxy = randn(1,2); # random walk
                Dxy = Dxy ./ norm(Dxy)
                V[id] = PC_v / 3 # slow speed
            end
        end
        # slow speed
        if ! isnan(Dmin) && Dmin<PC_h
            V[id] = PC_v / 5
        end
        # Store
        DXY[id,:] = Dxy
        DMIN[id] = Dmin
    end

end


#### HARVEST for a season
function fnc_harvest(school,fish,cons,fishtree,EVENTS,FLAGS);
    empty!(EVENTS["captor"]) #In this function we manage capture events
    empty!(EVENTS["captured"])
    
    CH,FX,TGT,DMIN=cons.H,fish.fx,cons.target,cons.Dmin
    
    nbcaptors=zeros(PC_n)
    Dharvest=PC_h #max harvest distance
    if FLAGS["implicit_fish"]
        Dhvarest=PF_sig+PC_h
    end
    
    Dmin=cons.Dmin
    f1=cons.measure["f1"]
    #For all fishers who have a target
    for i = EVENTS["targeting"]
        if Dmin[i]<Dharvest #If the target is within harvesting distance
            tgt=TGT[i]
            f1[i]+=1
            if FLAGS["implicit_fish"]
                #If the school is a simple disk with a population variable
                cons.V[i]=0
                school.pop[tgt]-=PC_q #use probability of catch as depletion rate
                cons.H[i]+=PC_q
                continue
            end
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
    for i=EVENTS["captor"]
        CH[i] += 1./nbcaptors[i]; #Fisher gets proportion of fish
    end
end



#### MOVE
#! CL <- school.x (school location); FX <- fish location;
#! distance
#! Randomly move the fish schools. 
#! Return the locations of the fish, school locations,
#! updated fisher locations
function fnc_move(school,fish,cons,fishtree,EVENTS,FLAGS)
    #CL=school.x;FX=fish.fx;FS=fish.fs;CC=cons.x;Dm=cons.Dmin;DXY=cons.DXY;
    CL,FX,CC,Dm,DXY,V=school.x,fish.fx,cons.x,cons.Dmin,cons.DXY,cons.V
    # schools and fish move
    empty!(EVENTS["jumped"])
    for i = 1:PS_n
        if school.pop[i]<1 || rand() < PS_p # if no fish left in school or randomly, jump
            j = school.fish[:,i];
            CL_x = rand() * GRD_mx; #center location
            CL_y = rand() * GRD_mx;
            CL[i,:] = [CL_x CL_y];
            school.pop[i]=length(j) #restore population counter to base level
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

    # fishers move
    # slow down as you approach fish
    #v = Dm; v[v.<PC_h] = PC_h; v[v.>PC_f] = PC_f;
    #v = (v .- PC_h) ./ (PC_f-PC_h);
    #range2 = PC_v - PC_vmn;
    #V = (v*range2) .+ PC_vmn;

    for f=1:PC_n
        CC_x = mod(CC[f,1] + (DXY[f,1]*V[f]), GRD_mx);
         CC_y = mod(CC[f,2] + (DXY[f,2]*V[f]), GRD_mx);
         cons.x[f,:] =[CC_x CC_y]
     end
     
     for f=EVENTS["targeting"]
         #reduce distance to target
         cons.Dmin[f]-=cons.V[f]
     end
    #return FX,CL,CC
end


