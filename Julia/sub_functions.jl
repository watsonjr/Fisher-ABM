
######## Functions for ABM
#### CONTACT NETWORK from social network
#! Iterates through social network adjacency matrix, generates 2 random numbers.
#! If random numbers are less than adjacency measure (friendship) they both
#! contact; else both do not contact. Return the network 2D array.
function fnc_contact(SN,CN)
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


#### DISTANCES
#!calculate distances between fish and fishermen? + search/steam switch
#! Returns Dx, the x-distances between fish and fishermen; Dy, y-distances;
#! D the Euclidean distances, MI which is a indicator for search steam switch (???)
function fnc_distance(FF,CC,MI)
    ## x,y, distances
    Dx = FF[:,1]' .- CC[:,1];
    Dy = FF[:,2]' .- CC[:,2];

    ## periodic boundaries
	GRD_mx_half = GRD_mx/2
    i = (abs(Dx).>(GRD_mx_half)) + (abs(Dy).>(GRD_mx_half));
    Dx[i.>1] = -sign(Dx[i.>1]).*(GRD_mx .- abs(Dx[i.>1]));
    Dy[i.>1] = -sign(Dy[i.>1]).*(GRD_mx .- abs(Dy[i.>1]));

    ## Distance
	@devec D = sqrt(Dx.^2 .+ Dy.^2); #! This line takes 40% of computation

	## search/steam switch
	MI[rand(PC_n).<=PC_rp] .-= 1;
	MI = abs(MI);

    return D,Dx,Dy,MI
end


#### INFORMATION and fisher direction
#! Iterate through the people you share information with and get the locations
#! of the fish within their view. Calculate the unit vector to the nearest fish
#! and if there's a fish in view, do a probabilistic catch. Otherwise roam around
#! randomly according to a self-correlated walk that approximates search behavior according to a Levy walk
#! Return the minimum distance, updated heading, JJ index of nearest fish
#! (KK) whether or not you caught something (?)

function fnc_information(D,Dx,Dy,DXY,MI,CN,id)
    # Who you will get info from
    # get the vector of people with whom you are currently in contact
    j = find(CN[id,:].==1); #!OPT devectorize this 

    # get information from friends and yourself
    DDi = 0;
    for k = 1:length(j) #!through friends list
		d = D[j[k],:]; #!return vector of distances for contact k
		kk = find(d.<PC_f); # find only those in fish finder view
		DDi = hcat(DDi,kk');
    end
    DDi = DDi[2:end];
    DD  = D[id,DDi];

    # find nearest fish
    if isempty(DD)==0 # if I see something
        Dmin = minimum(DD);
        II = find(DD .== Dmin); #!find index of lowest distance?
        JJ = DDi[II][1];

        # calculate unit vector DXY to nearest fish
        DXY = [Dx[id,JJ][1] Dy[id,JJ][1]];
		@devec DXY = DXY ./ sqrt(DXY[1].^2 + DXY[2].^2);

        if Dmin <= PC_h #if there's a fish within view
            # probabilistic catch; if successful, then catch fish
            r =rand();if r < PC_q;KK = 1;else;KK = 0;end
        else
            KK = 0; #did you harvest or not?
        end

    else # else I roam around randomly 
        Dmin = 999; # flag, if you don't see anything
        JJ = 0; KK = 0;

        #Steam or search pattern
        if MI == 0 #
			DXY = DXY + (randn(1,2).*PC_r1)
        	DXY = DXY ./ sqrt(DXY[1].^2 + DXY[2].^2);
		else
			DXY = DXY + (randn(1,2).*PC_r2)
        	DXY = DXY ./ sqrt(DXY[1].^2 + DXY[2].^2);
        end
    end

    return Dmin,DXY,JJ,KK
end


#### HARVEST for a season
function fnc_harvest(KK,JJ,CC,CS,FF);
    II    = KK.*JJ
    IIu   = unique(II);

    for i in IIu
        j = find(II.==i);
        CC[j] =  (1 / length(j)) * KK[j];
    end
    FF[IIu[IIu.>0],:] = NaN;
    CS = CS + CC;
    return CC,CS,FF
end


#### MOVE
#! 6/11 implementation: CL <- fish.sx (school location); FX <- fish location;
#! FS <- index of the fish's school; CC <- fisher location?; D1 <- cumulative distance
#! Randomly move the fish schools. Return the locations of the fish, school locations,
#! updated fisher locations
function fnc_move(CL,FX,FS,CC,DXY,D1)

	# schools move randomly
    rn = rand(PS_n);
    i = find(rn .<= PS_p);
    CL_x = CL[:,1]; CL_y = CL[:,2];
    CL_x[i] = rand() * GRD_mx;
    CL_y[i] = rand() * GRD_mx;

    # fish move
    for j = 1:length(i) #! for each fish in the school that moves
        Fx = mod(CL_x[i[j]].+(randn(PF_n,1)*PF_sig),GRD_mx); #! move fish with periodic boundary
    	Fy = mod(CL_y[i[j]].+(randn(PF_n,1)*PF_sig),GRD_mx);
		FX[FS.==i[j],1] = Fx; #! update fish locations
		FX[FS.==i[j],2] = Fy;
	end

    # fishers move
    CC_x = mod(CC[:,1] .+ (DXY[:,1].*PC_v), GRD_mx);
    CC_y = mod(CC[:,2] .+ (DXY[:,2].*PC_v), GRD_mx);
    D2 = D1.+1; #update cumulative distance

    return FX,[CL_x CL_y],[CC_x CC_y], D2
end

### Connectance
function fnc_connectance(SN)
	C = sum(SN,1);
	C = sum(C,2);
	C = squeeze(C,1);
	C = squeeze(C,1);
	return C
end


