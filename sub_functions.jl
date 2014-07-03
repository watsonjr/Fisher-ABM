
######## Functions for ABM
#### CONTINUOUS DISTANCES
#!calculate distances between fish and fishermen? + search/steam switch
#! Returns Dx, the x-distances between fish and fishermen; Dy, y-distances;
#! D the Euclidean distances, MI which is a indicator for search steam switch (???)
#function fnc_distance(FF,CC,MI)
#    ## x,y, distances
#    Dx = FF[:,1]' .- CC[:,1];
#    Dy = FF[:,2]' .- CC[:,2];
#
#    ## periodic boundaries
#	GRD_mx_half = GRD_mx/2
#    i = (abs(Dx).>(GRD_mx_half)) + (abs(Dy).>(GRD_mx_half));
#    Dx[i.>1] = -sign(Dx[i.>1]).*(GRD_mx .- abs(Dx[i.>1]));
#    Dy[i.>1] = -sign(Dy[i.>1]).*(GRD_mx .- abs(Dy[i.>1]));
#
#    ## Distance
#	@devec D = sqrt(Dx.^2 .+ Dy.^2); #! This line takes 40% of computation
#
#	## search/steam switch
#	MI[rand(PC_n).<=PC_rp] .-= 1;
#	MI = abs(MI);
#
#    return D,Dx,Dy,MI
#end

#### DISCRETE DISTANCES
#!calculate distances between fish and fishermen? + search/steam switch
#! Returns Dx, the x-distances between fish and fishermen; Dy, y-distances;
#! D the Euclidean distances, MI which is a indicator for search steam switch (???)
function fnc_fishfinder(Fx,Cx,grd,PC_f)
	Fi = round(Fx) .+ 1;
	Ci = round(Cx) .+ 1;
	Ni = fill(0,PC_n,1);

	for i = 1:PC_n
		# get index of all fish near fisher i
		DX = Ci[i,1] .- Fi[:,1]
		DY = Ci[i,2] .- Fi[:,2]
		j = (abs(DX).>grd) + (abs(DY).>grd);
		DX[j.>1] = -sign(DY[j.>1]).*(GRD_mx .- abs(DX[j.>1]));
		DY[j.>1] = -sign(DX[j.>1]).*(GRD_mx .- abs(DY[j.>1]));
		I  = find(((abs(DX) .< PC_f) + (abs(DY) .< PC_f)) .== 2);

		# find squared euclidean distance with these nearish fish
		if isempty(I)==0
			dx = Fx[I,1] .- Cx[i,1];
			dy = Fx[I,2] .- Cx[i,2];
			j = (abs(dx).>grd) + (abs(dy).>grd);
			dx[j.>1] = -sign(dx[j.>1]).*(GRD_mx .- abs(dx[j.>1]));
			dy[j.>1] = -sign(dy[j.>1]).*(GRD_mx .- abs(dy[j.>1]));
			D   = abs2(dx) + abs2(dy);
		 	dmi = find(D.==minimum(D));

			# store nearest fish
			Ni[i] = int(I[dmi[1]]); # Index of nearest fish in fishfinder
		end
	end

    return Ni # Index of nearest fish for each fisher (in fishfinder)
end


##### Distance calculation using KDTREE
#function fnc_distance_3(FF,CC,Fish_finder)
#	
#	cons_Ni = Array(Int,PC_n,1)
#	cons_Nd = Array(Float64,PC_n,1)
#	Dxy     = Array(Float64,PC_n,2)
#	for i = 1:PC_n
#		ni,nd = nearest(FTREE,squeeze(CC[i,:]',2),1);
#		if nd[1] > Fish_finder; nd[1]=NaN; ni[1]=0; end
#		cons_Ni[i] = ni[1];
#		cons_Nd[i] = nd[1];
#	end
#
#    return cons_Ni, cons_Nd
#end

#### Search/steam switch
#! that is, is a fisher can't see any fish
#! they can either spin around or move in a straightish line
#! they switch between this probabilistically
function fnc_steam(MI)
	MI[rand(PC_n).<=PC_rp] .-= 1;
	MI = abs(MI);
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
	dxy=cons.DXY;Ni=cons.Ni;Fx=fish.fx;Cx=cons.x;MI=cons.MI;

	DXY  = Array(Float64,PC_n,2) # heading
	DMIN = Array(Float64,PC_n) # shortest distance
	JJ   = fill(0,PC_n) # index of nearest fish
	KK   = fill(0,PC_n) # index of whether you caught fish

	for id = 1:PC_n
		# get the vector of people with whom you are currently in contact
		J  = find(CN[id,:].==1); # index of friends
		Jn = length(J); # number of friends 

		# calculate distances to all fish you have info on
		DD = fill(NaN,PC_n) # distance to your fish and friend's
		dx = fill(NaN,PC_n,2) # direction to fish
		for i = 1:Jn # for each friend
			ii = J[i]; # get his/her index
			if Ni[ii] != 0 # if they see a fish
				dx[i,:] = Fx[Ni[ii],:] - Cx[id,:]; # calc you direction to it
				DD[i] = sqrt(abs2(dx[i,1]) + abs2(dx[i,2])); # and distance
			end
		end
		Dmin = minimum(DD); # shortest distance to a fish

		# Action-decide heading, catch fish
		if isnan(Dmin) == 0 # if I see anything
			ii = find(DD .== Dmin); #!index of friend next to fish
			jj = Ni[J[ii]]; #!index of nearest fish 

			# calculate unit vector DXY to nearest fish
			Dxy = dx[ii,:] ./ norm(dx[ii,:]); 

			if Dmin <= PC_h #if there's a fish within view
				# probabilistic catch; if successful, then catch fish
				r =rand();if r < PC_q;kk = 1;else;kk = 0;end
			else
				kk = 0; #index of whther harvest or not?
			end

		else # else I roam around randomly 
			Dmin = 999; # flag, if you don't see anything
			jj = 0; kk = 0;

			#Steam or search pattern
			if MI == 0 #
				Dxy = dxy[id,:] + (randn(1,2).*PC_r1)
				Dxy = Dxy ./ norm(Dxy)
			else
				Dxy = dxy[id,:] + (randn(1,2).*PC_r2)
				Dxy = Dxy ./ norm(Dxy)
			end
		end

		# Store
		DXY[id,:] = Dxy
		DMIN[id] = Dmin
		JJ[id] = jj[1];
		KK[id] = kk[1];
	end

    return DMIN,DXY,JJ,KK
end


#### HARVEST for a season
#! KK is the index of whether a given fish has caught a fish
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
function fnc_move(CL,FX,FS,CC,Dm,DXY)

	# schools move randomly
    rn = rand(PS_n);
    i = find(rn .<= PS_p);
    CL_x = CL[:,1]; CL_y = CL[:,2];
    CL_x[i] = rand() * GRD_mx;
    CL_y[i] = rand() * GRD_mx;

    # fish move
    for j = 1:length(i) #! for each fish in the school that moves
        Fx = mod(CL_x[i[j]].+(randn(PF_n,1)*PF_sig),GRD_mx);
    	Fy = mod(CL_y[i[j]].+(randn(PF_n,1)*PF_sig),GRD_mx);
		FX[FS.==i[j],1] = Fx; #! update fish locations
		FX[FS.==i[j],2] = Fy;
	end

    # fishers move
    v = Dm; v[v.<2] = 2; v[v.>PC_f] = PC_f; # slow down as you approach fish
    v = (v .- 2) ./ (10-2);
    range2 = PC_v - 0.2;
	V = (v*range2) .+ 0.2;

    CC_x = mod(CC[:,1] .+ (DXY[:,1].*V), GRD_mx);
    CC_y = mod(CC[:,2] .+ (DXY[:,2].*V), GRD_mx);

    return FX,[CL_x CL_y],[CC_x CC_y]
end


