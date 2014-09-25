
######## Functions for ABM

#### HAUL TIME
#! running time between hauls
#! and estimate the running mean time between schools
#! and estimate the difference in this running mean 
#! which is the switch for the while loop
function fnc_tau(dTs,cons,FLAGS)
	Ts,Tv,ts,ns=cons.Ts,cons.Tv,cons.ts,cons.ns
	#H=cons.H;Ts=cons.Ts;ts=cons.ts;sn=cons.sn;
	for I = 1:PC_n
		# if you caught fish
		# and the cumulative haul time (ts) is large 
		# then you've encountered a new school
		if in(I,FLAGS["captor"])
			if ts[I] > (10*PF_sig/PC_v) 
				ns[I] += 1; # update school counter
				Ts_old = Ts[I]; # current mean
				Tv_old = Tv[I]; # current variance

				Ts[I] = Ts_old + ((ts[I]-Ts_old)/ns[I]) # run mean
				Tv[I] = Tv_old + ((ts[I]-Ts_old)*(ts[I]-Ts[I])) # run variance
				
				dTs[I] = abs(Ts[I]-Ts_old)/Ts[I]; # fractional change in mean
				ts[I] = 1; # reset how long it took to find school
			else
				ts[I] = 1;
			end
		else # if you didn't catch anythin
			ts[I] += 1
		end
	end
	#return Ts,Tv,ts,ns,dTs
end


#### FISH FINDER / DISTANCE
#!calculate distances between fish and fishermen? + search/steam switch
#! Returns Dx, the x-distances between fish and fishermen; Dy, y-distances;
#! D the Euclidean distances, MI which is a indicator for search steam switch (???)
function fnc_fishfinder(grd,PC_f,school,fish,cons,FLAGS)
	Fx,Sx,Si,Cx,fish_per_school=fish.fx,school.x,fish.fs,cons.x,school.fish
	# First, find fishers that are likely near fish
	# by gauging distance to all school centres
	II = cell(PC_n); # index of schools each fisher is near
	for i = 1:PC_n
		II[i] =[]
		if in(true,FLAGS["benichou"]) &&   cons.MI[i]!=1 #! ( cons.MI[i]=1 ||cons.Ni[i]==0 || isnan(fish.fx[cons.Ni[i]]) )
			continue
		end
		for j=1:PS_n
			D=hypot( (Cx[i,:]-Sx[j,:] )...)
			if D.<((2.5.*PF_sig).+PC_f)
				II[i]=[II[i],j]
			end
		end
	end

	empty!(FLAGS["new_neighbor"])
	# Output - index of nearest fish for all fishers (if any)
	Ni = cons.Ni;
	# find for each fisher
	for i = 1:PC_n

		if isempty(II[i]) == 0 # if there is a school nearby fisher i

			# if so, get index k of fish in all schools, near fisher i
			newk=0
			Dmin=NaN #minimal distance to another fish
			for j = 1:length(II[i])
				for k = fish_per_school[:,II[i][j]]
					D=hypot( (Cx[i,:]-Fx[k,:] )%(GRD_mx2)...)
					if D<PC_f && (isnan(Dmin) || D<Dmin)
						newk=k
						Dmin=D
					end
				end
			end
			if newk!=Ni[i]
				Ni[i]=newk
				if newk!=0
					push!(FLAGS["new_neighbor"],i)
				end
			end
		else 
				Ni[i] = 0;
		end
	end

	#cons.Ni=Ni # Index of nearest fish for each fisher (in fishfinder)
	#return Ni
end


#### spatial Difference function
#! x1 = first x,y location
#! x2 = second x,y location
#! dx,dy = difference in x,y accounting for periodic boundary
function fnc_difference(x1,x2)
	# difference
	dx = x1[:,1] .- x2[:,1]
	dy = x1[:,2] .- x2[:,2]
	# periodic boundary
	j = (abs(dy).>GRD_mx2) + (abs(dx).>GRD_mx2);
	dx[j.>1] = -sign(dx[j.>1]).*(GRD_mx .- abs(dx[j.>1]));
	dy[j.>1] = -sign(dy[j.>1]).*(GRD_mx .- abs(dy[j.>1]));

	return dx,dy
end



#### Search/steam switch
#! that is, is a fisher can't see any fish
#! they can either spin around or move in a straightish line
#! they switch between this probabilistically
function fnc_steam(school,fish,cons,FLAGS)
	MI=cons.MI
	for i = 1:PC_n
		if MI[i] == 1 # if steaming
			if rand() < PC_rp # maybe switch to tumbling
				MI[i] = 0
			end
		else # if tumbling
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
function fnc_contact(school,fish,cons,FLAGS)
	SN=cons.SN
	CN = zeros(Int,PC_n,PC_n)
	probing= FLAGS["new_neighbor"]  #ask to any friend who has a new neighbor
	#if !isempty(probing)
	#	print(probing)
	#end
	for i = 1:PC_n
		#for j = i:PC_n
		for j in probing
			f1 = SN[i,j];
			f2 = SN[j,i];
			RN = rand(2,1);

			if RN[1] < f1 && RN[2] < f2
				CN[i,j] = 1;
				CN[j,i] = 1;
			else
				CN[i,j] = 0;
				CN[j,i] = 0;
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
function fnc_decision(source,target,school,fish,cons,FLAGS)
	return 1
end

#### INFORMATION and fisher direction
#! Iterate through the people you share information with and get the locations
#! of the fish within their view. Calculate the unit vector to the nearest fish.
#! Otherwise roam around randomly according to a self-correlated walk that 
#! approximates search behavior according to a Levy walk
#! Return the minimum distance, updated heading, JJ index of nearest fish
function fnc_information(CN,school,fish,cons,FLAGS)
	#In this function we manage targeting events
	
	dxy,Ni,Fx,Cx,MI=cons.DXY,cons.Ni,fish.fx,cons.x,cons.MI

	DXY  = Array(Float64,PC_n,2) # heading
	DMIN = Array(Float64,PC_n) # shortest distance
	V	 = cons.V # speed

	#Fishers for whom previous target has been canceled:
	#fish captured or part of school that jumped
	for id=FLAGS["targeting"]
		tgt=cons.target[id]
		if in(tgt, FLAGS["captured"]) || in(fish.fs[tgt], FLAGS["jumped"]) 
			cons.Dmin[id]=NaN
			cons.target[id]=0
			delete!(FLAGS["targeting"],id)
			V[id]=PC_v
		end
	end
	
	for id = 1:PC_n
		#if in(true,FLAGS["benichou"])& MI[id] == 0
		#	continue
		#end
		# get the vector of people with whom you are currently in contact
		J  = find(CN[id,:].==1); # index of friends
		Jn = length(J); # number of friends 

		# calculate distances to all fish you have info on
		DD = fill(NaN,PC_n) # distance to your fish and friend's
		Dx = fill(NaN,PC_n) # dx
		Dy = fill(NaN,PC_n) # dy
		for i = 1:Jn # for each friend
			ii = J[i]; # get his/her index
			if Ni[ii] != 0 && !isnan(Fx[Ni[ii],1]) # if they see a fish
				(dx,dy) = fnc_difference(Fx[Ni[ii],:],Cx[id,:]);
				DD[i] = sqrt(abs2(dx) + abs2(dy))[1]; # and distance
				Dx[i] = dx[1]; Dy[i] = dy[1];
			end
		end
		Dmin = minimum(DD); # shortest distance to a fish
		current=cons.Dmin[id]
		
		# Action-decide heading, catch fish
		if  !isnan(Dmin)  && (isnan(current) || Dmin<current ) 
			# if I see anything better than current target
			ii = find(DD .== Dmin); #!index of friend next to fish
			ii = ii[1]; # if more than one friend is next to the same fish
			jj = Ni[J[ii]]; #!index of nearest fish 

			# calculate unit vector DXY to nearest fish
			Dxy = [Dx[ii] Dy[ii]] ./ norm([Dx[ii] Dy[ii]]); 

			push!(FLAGS["targeting"],id)
			cons.target[id]=jj
				
		elseif  in(id,FLAGS["targeting"]) && !isnan(current) #Stay with current target
			#print( id," ",current," ", Dmin,"\n")
			Dxy=cons.DXY[id,:]
			Dmin=cons.Dmin[id]
		else # else I roam around randomly 
			cons.target[id]=0
			delete!(FLAGS["targeting"],id)

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

	cons.Dmin,cons.DXY= DMIN,DXY
end


#### HARVEST for a season
function fnc_harvest(school,fish,cons,FLAGS);
	empty!(FLAGS["captor"]) #In this function we manage capture events
	empty!(FLAGS["captured"])
	
	CH,FX,TGT,DMIN=cons.H,fish.fx,cons.target,cons.Dmin
	
	#For all fishers that have a target
	nbcaptors=zeros(PC_n)
	for i = FLAGS["targeting"]
		if cons.Dmin[i]<PC_h #If the target is within harvesting distance
			cons.V[i]=0 #The boat stops moving until it catches the fish
			# probabilistic catch
			r =rand()
			if r < PC_q
				push!(FLAGS["captured"],TGT[i])
				FX[TGT[i],:]=NaN; #Fish disappears
				push!(FLAGS["captor"],i)
				nbcaptors[i]+=1
				cons.V[i] = PC_v
			end
		end
	end
	for i=FLAGS["captor"]
		CH[i] += 1./nbcaptors[i]; #Fisher gets proportion of fish
	end
end
#= Previous version: 
#! KK is the index of whether a given fisher has caught a fish
#! JJ is the index of the fish he's caught
#! CH is the cumulative catch
#! FF are the fish locations, which must be updated is fish are caught
	CH,FF=cons.H,fish.fx
	II    = KK.*JJ
	IIu   = unique(II);

	for i in IIu
		j = find(II.==i);
		CH[j] = CH[j] + ((1 / length(j)) * KK[j]);
	end
	FF[IIu[IIu.>0],:] = NaN;
	#return CH,FF
end =#


#### MOVE
#! CL <- school.x (school location); FX <- fish location;
#! distance
#! Randomly move the fish schools. 
#! Return the locations of the fish, school locations,
#! updated fisher locations
function fnc_move(school,fish,cons,FLAGS)
	#CL=school.x;FX=fish.fx;FS=fish.fs;CC=cons.x;Dm=cons.Dmin;DXY=cons.DXY;
	CL,FX,CC,Dm,DXY,V=school.x,fish.fx,cons.x,cons.Dmin,cons.DXY,cons.V
	# schools and fish move
	empty!(FLAGS["jumped"])
	for i = 1:PS_n
		j = school.fish[:,i];
		k = FX[j,1];
		m = find(isnan(k).==0)
		if isempty(m) == 1 || rand() < PS_p # if no fish left in school or randomly, jump
			CL_x = rand() * GRD_mx;
			CL_y = rand() * GRD_mx;
			F_x  = mod(CL_x.+(randn(PF_n,1)*PF_sig),GRD_mx);
			F_y  = mod(CL_y.+(randn(PF_n,1)*PF_sig),GRD_mx);
			CL[i,:] = [CL_x CL_y];
			FX[j,:] = [F_x F_y];
			push!(FLAGS["jumped"],i);
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
 	
 	for f=FLAGS["targeting"]
 		#reduce distance to target fish
 		cons.Dmin[f]-=cons.V[f]
 	end
    #return FX,CL,CC
end


