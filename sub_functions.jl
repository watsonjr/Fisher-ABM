
######## Functions for ABM

#### HAUL TIME
#! running time between hauls
#! and estimate the running mean time between schools
#! and estimate the difference in this running mean 
#! which is the switch for the while loop
function fnc_tau(H,Ts,Tv,ts,ns,dTs)
	#H=cons.H;Ts=cons.Ts;ts=cons.ts;sn=cons.sn;
	for I = 1:PC_n
		# if you caught fish
		# and the cumulative haul time (ts) is large 
		# then you've encountered a new school
		if H[I] == 1
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
	return Ts,Tv,ts,ns,dTs
end


#### FISH FINDER / DISTANCE
#!calculate distances between fish and fishermen? + search/steam switch
#! Returns Dx, the x-distances between fish and fishermen; Dy, y-distances;
#! D the Euclidean distances, MI which is a indicator for search steam switch (???)
function fnc_fishfinder(Fx,Sx,Si,Cx,grd,PC_f)

	# First, find fishers that are likely near fish
	# by gauging distance to all school centres
	II = cell(PC_n); # index of schools each fisher is near
	for i = 1:PC_n
		(dx,dy) = fnc_difference(Cx[i,:],Sx);
		D = sqrt(abs2(dx) .+ abs2(dy));
		JJ = find(D.<((2.5.*PF_sig).+PC_f));
		II[i] = JJ # index of nearby schools (empty if none)
	end

	# Output - index of nearest fish for all fishers (if any)
	Ni = fill(0,PC_n,1);

	# find for each fisher
	for i = 1:PC_n

		if isempty(II[i]) == 0 # if there is a school nearby fisher i

			# if so, get index k of fish in all schools, near fisher i
			k = Array(Float64,0)
			for j = 1:length(II[i])
				k = [k, find(Si .== II[i][j])];
			end

			# get dx and dy,
			Ci = round(Cx[i,:]) .+ 1;
			Fi = round(Fx[k,:]) .+ 1;
			(dx,dy) = fnc_difference(Ci,Fi)

			# get index of those fish, in schools k, near fisher i
			I  = find(((abs(dx) .< PC_f) + (abs(dy) .< PC_f)) .== 2);

			# find squared euclidean distance with these nearish fish
			if isempty(I)==0
				(dx,dy) = fnc_difference(Fx[k[I],:],Cx[i,:]);
				D   = abs2(dx) + abs2(dy);
				dmi = find(D.==minimum(D));

				# store nearest fish
				Ni[i] = int(k[I[dmi[1]]]); # Index of nearest fish in fishfinder
			end

		else 
				Ni[i] = 0;
		end
	end

    return Ni # Index of nearest fish for each fisher (in fishfinder)
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
function fnc_steam(MI)
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
function fnc_contact(SN)
	CN = Array(Int,PC_n,PC_n)
	for i = 1:PC_n
		for j = i:PC_n
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


#### INFORMATION and fisher direction
#! Iterate through the people you share information with and get the locations
#! of the fish within their view. Calculate the unit vector to the nearest fish
#! and if there's a fish in view, do a probabilistic catch. Otherwise roam around
#! randomly according to a self-correlated walk that approximates search behavior according to a Levy walk
#! Return the minimum distance, updated heading, JJ index of nearest fish
#! (KK) whether or not you caught something (?)
function fnc_information(dxy,Ni,Fx,Cx,MI,CN)
	#dxy=cons.DXY;Ni=cons.Ni;Fx=fish.fx;Cx=cons.x;MI=cons.MI;

	DXY  = Array(Float64,PC_n,2) # heading
	DMIN = Array(Float64,PC_n) # shortest distance
	V	 = Array(Float64,PC_n) # speed
	JJ   = fill(0,PC_n) # index of nearest fish
	KK   = fill(0,PC_n) # index of whether you caught fish

	for id = 1:PC_n
		# get the vector of people with whom you are currently in contact
		J  = find(CN[id,:].==1); # index of friends
		Jn = length(J); # number of friends 

		# calculate distances to all fish you have info on
		DD = fill(NaN,PC_n) # distance to your fish and friend's
		Dx = fill(NaN,PC_n) # dx
		Dy = fill(NaN,PC_n) # dy
		for i = 1:Jn # for each friend
			ii = J[i]; # get his/her index
			if Ni[ii] != 0 # if they see a fish
				(dx,dy) = fnc_difference(Fx[Ni[ii],:],Cx[id,:]);
				DD[i] = sqrt(abs2(dx) + abs2(dy))[1]; # and distance
				Dx[i] = dx[1]; Dy[i] = dy[1];
			end
		end
		Dmin = minimum(DD); # shortest distance to a fish

		# Action-decide heading, catch fish
		if isnan(Dmin) == 0 # if I see anything
			ii = find(DD .== Dmin); #!index of friend next to fish
			ii = ii[1]; # if more than one friend is next to the same fish
			jj = Ni[J[ii]]; #!index of nearest fish 

			# calculate unit vector DXY to nearest fish
			Dxy = [Dx[ii] Dy[ii]] ./ norm([Dx[ii] Dy[ii]]); 

			if Dmin <= PC_h #if there's a fish within view
				# probabilistic catch; if successful, then catch fish
				r =rand();if r < PC_q;kk = 1;else;kk = 0;end
			else
				kk = 0; #index of whther harvest or not?
			end

			# slow speed
			V[id] = PC_v / 5
				
		else # else I roam around randomly 
			Dmin = 999; # flag, if you don't see anything
			jj = 0; kk = 0;

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

		# Store
		DXY[id,:] = Dxy
		DMIN[id] = Dmin
		JJ[id] = jj[1];
		KK[id] = kk[1];
	end

    return DMIN,DXY,JJ,KK,V
end


#### HARVEST for a season
#! KK is the index of whether a given fisher has caught a fish
#! JJ is the index of the fish he's caught
#! CH is the cumulative catch
#! FF are the fish locations, which must be updated is fish are caught
function fnc_harvest(KK,JJ,CH,FF);
    II    = KK.*JJ
    IIu   = unique(II);

    for i in IIu
        j = find(II.==i);
        CH[j] = CH[j] + ((1 / length(j)) * KK[j]);
    end
    FF[IIu[IIu.>0],:] = NaN;
    return CH,FF
end


#### MOVE
#! CL <- fish.sx (school location); FX <- fish location;
#! distance
#! Randomly move the fish schools. 
#! Return the locations of the fish, school locations,
#! updated fisher locations
function fnc_move(CL,FX,FS,CC,Dm,DXY,V)
	#CL=fish.sx;FX=fish.fx;FS=fish.fs;CC=cons.x;Dm=cons.Dmin;DXY=cons.DXY;

	# schools and fish move
	for i = 1:PS_n
		j = find(FS.==i);
		k = FX[j,1];
		m = find(isnan(k).==0)
		if isempty(m) == 1# if no fish left in school
			CL_x = rand() * GRD_mx;
			CL_y = rand() * GRD_mx;
			F_x  = mod(CL_x.+(randn(PF_n,1)*PF_sig),GRD_mx);
			F_y  = mod(CL_y.+(randn(PF_n,1)*PF_sig),GRD_mx);
			CL[i,:] = [CL_x CL_y];
			FX[j,:] = [F_x F_y];
		elseif rand() < PS_p # else maybe jump
			CL_x = rand() * GRD_mx;
			CL_y = rand() * GRD_mx;
			F_x  = mod(CL_x.+(randn(PF_n,1)*PF_sig),GRD_mx);
			F_y  = mod(CL_y.+(randn(PF_n,1)*PF_sig),GRD_mx);
			CL[i,:] = [CL_x CL_y];
			FX[j,:] = [F_x F_y];
		end
	end

    # fishers move
    # slow down as you approach fish
    #v = Dm; v[v.<PC_h] = PC_h; v[v.>PC_f] = PC_f;
    #v = (v .- PC_h) ./ (PC_f-PC_h);
    #range2 = PC_v - PC_vmn;
	#V = (v*range2) .+ PC_vmn;

    CC_x = mod(CC[:,1] .+ (DXY[:,1].*V), GRD_mx);
    CC_y = mod(CC[:,2] .+ (DXY[:,2].*V), GRD_mx);
    CC = [CC_x CC_y];

    return FX,CL,CC
end


