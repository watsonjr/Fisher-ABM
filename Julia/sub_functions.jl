
######## Functions for ABM
#### CONTACT NETWORK from social network
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
function fnc_distance(FF,CC,MI)
    ## x,y, distances
    Dx = FF[:,1]' .- CC[:,1];
    Dy = FF[:,2]' .- CC[:,2];

    ## periodic boundaries
    i = (abs(Dx).>(GRD_mx/2)) + (abs(Dy).>(GRD_mx/2));
    Dx[i.>1] = -sign(Dx[i.>1]).*(GRD_mx .- abs(Dx[i.>1]));
    Dy[i.>1] = -sign(Dy[i.>1]).*(GRD_mx .- abs(Dy[i.>1]));

    ## Distance
    D = sqrt(Dx.^2 + Dy.^2);

	# search/stream switch
	MI[rand(PC_n).<=PC_rp] .-= 1;
	MI = abs(MI);

    return D,Dx,Dy,MI
end


#### INFORMATION and fisher direction
function fnc_information(D,Dx,Dy,DXY,MI,CN,id)
    # Who you will get info from
    j = find(CN[id,:].==1);

    # get information from friends and yourself
    DDi = 0;
    for k = 1:length(j)
        d = D[j[k],:];
        kk = find(d.<PC_f); # find only those in fish finder view
        DDi = hcat(DDi,kk');
    end
    DDi = DDi[2:end];
    DD  = D[id,DDi];

    # find nearest fish
    if isempty(DD)==0 # if I see something
        Dmin = minimum(DD);
        II = find(DD .== Dmin);
        JJ = DDi[II][1];
        DXY = [Dx[id,JJ][1] Dy[id,JJ][1]];
        DXY = DXY ./ sqrt(DXY[1].^2 + DXY[2].^2);

        if Dmin <= PC_h
            # probabilistic catch; if successful, then catch fish
            r =rand();if r < PC_q;KK = 1;else;KK = 0;end
        else
            KK = 0;
        end

    else # else I roam around randomly 
        Dmin = 999; # flag, if you don't see anything
        JJ = 0; KK = 0;

       	if MI == 0
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
function fnc_move(CL,FX,FS,CC,DXY,D1)

    # schools move randomly
    rn = rand(PS_n);
    i = find(rn .<= PS_p);
    CL_x = CL[:,1]; CL_y = CL[:,2];
    CL_x[i] = rand() * GRD_mx;
    CL_y[i] = rand() * GRD_mx;

    # fish move
    for j = 1:length(i)
    	Fx = mod(CL_x[i[j]].+(randn(PF_n,1)*PF_sig),GRD_mx);
    	Fy = mod(CL_y[i[j]].+(randn(PF_n,1)*PF_sig),GRD_mx);
		FX[FS.==i[j],1] = Fx;
		FX[FS.==i[j],2] = Fy;
	end

    # fishers move
    CC_x = mod(CC[:,1] .+ (DXY[:,1].*PC_v), GRD_mx);
    CC_y = mod(CC[:,2] .+ (DXY[:,2].*PC_v), GRD_mx);
    D2 = D1.+1;

    return FX,[CL_x CL_y],[CC_x CC_y], D2
end

function connectance(SN)
	C = sum(SN,2);
	C = sum(C,1);
	C = squeeze(C,2);
	C = squeeze(C,1);
	return C
end

