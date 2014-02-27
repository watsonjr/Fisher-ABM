
######## Functions for ABM


#### FISH RELOCATION
function fnc_relocate(Cl,FF)
    i = int(ceil(rand(1)*(P_F_m.*P_F_n))); # choose a number of fish to relocate
    j = shuffle([1:P_F_n]); # shuffle for randomness
    j = j[1:i[1]]; # choose the select few
    k = ceil(rand(i[1],1) * P_Cl_n); # cluster center to move to
    a = repmat([Grd_mx Grd_my],length(j),1); # for periodic boundaries
    F = mod(Cl[k[:],:] + randn(i[1],2).*P_F_dx, a); # move them
    FF[j,:] = F; # add back
    return FF
end


#### DISTANCES
function fnc_distance(FF,CC)
    ## x,y, distances
    Dx = FF[:,1]' .- CC[:,1];
    Dy = FF[:,2]' .- CC[:,2];

    ## periodic boundaries
    i = (abs(Dx).>(Grd_mx/2)) + (abs(Dy).>(Grd_my/2));
    Dx[i.>1] = -sign(Dx[i.>1]).*(Grd_mx - abs(Dx[i.>1]));
    Dy[i.>1] = -sign(Dy[i.>1]).*(Grd_my - abs(Dy[i.>1]));

    ## Distance
    D = sqrt(Dx.^2 + Dy.^2);
    return D,Dx,Dy
end


#### INFORMATION
function fnc_information(D,Dx,Dy,id)
    # who is in my social network
    k = P_C_Sn[id,:];
    RN = rand(1,P_C_n); # bernoulli probabilistic
    j = find(RN.<k);

    # get information from friends and yourself
    DD = 0; DDx = 0; DDy = 0;
    for k = 1:length(j)
        d = D[j[k],:];
        kk = find(d.<P_C_ff); # find only those in fish finder view
        DD  = hcat(DD,D[j[k],kk]); ### Choose only those fish in finder
        DDx = hcat(DDx,Dx[j[k],kk]);
        DDy = hcat(DDy,Dy[j[k],kk]);
    end
    DD  = DD[2:end];
    DDx = DDx[2:end];
    DDy = DDy[2:end];

    # find nearest fish
    if isempty(DD)==0 # if I see something
        Dmin = min(DD);
        II = find(D .== Dmin);
        K,JJ = ind2sub(size(D),II); # index of nearest fish (JJ)
        JJ = JJ[1];
        DDx = Dx[id,JJ][1];
        DDy = Dy[id,JJ][1];
    else # else I roam around randomly
        Dmin = Grd_mx;
        JJ = 0;
        DDx  = 0; rand();
        DDy  = 0; rand();
    end

    return Dmin,DDx,DDy,JJ
end


#### DIRECTION and SPEED
function fnc_direction(Dmin,DDx,DDy,ANG_in)
    if Dmin > P_C_c && Dmin <= P_C_ff # if nearish move towards
        DXY  = [DDx DDy] ./ Dmin;
        ANG_out  = DXY[:,2] ./ DXY[:,1];
        vr   = (Dmin ./ P_C_ff) .* P_C_v;
        KK   = 0;
    elseif Dmin <= P_C_c # if very near - go fishing
        DXY = [DDx DDy] ./ Dmin;
        ANG_out  = DXY[:,2] ./ DXY[:,1];
        vr = 0; # no speed, stay where you are

        # probabilistic catch
        r =rand();
        if r < P_C_pr # if successful, then catch fish
            KK  = 1;
        else
            KK = 0;
        end
    else # if not so near, random walk
        ANG_out = ANG_in + (P_C_ang * rand() - (P_C_ang/2));
        vr   = 1.*P_C_v;
        KK   = 0;
    end
    return ANG_out[1],vr,KK
end


#### HARVEST
function fnc_harvest(KK,JJ,CC,FF,CLi,CLx)
    II = KK.*JJ
    if II != 0 # if the fish is caught
        # record that is it
        CC = 1; # catch
        # and relocate fish (to far away cl centre picked at random)
        LL = CLi[II];
        id = P_Cl_id[P_Cl_id.!=LL];
        id = shuffle(id);
        FF[II,:] = mod(CLx[id[1],:] + (randn().*P_F_dx),[Grd_mx Grd_my]);
    else
        CC = 0; # no catch
    end
    return CC, FF
end


#### MOVE
function fnc_move(CL,CC,theta_c,VR)

    # clusters centers move
    ## Make this probabilistic (i.e. doesn't move all the time)
    theta_f = rand(P_Cl_n)*2*pi; # angle
    f = rand(P_Cl_n,1).^(-1/P_Cl_al); # magnitude
    CL_x = mod(CL[:,1] + (f.*cos(theta_f)), Grd_mx);
    CL_y = mod(CL[:,2] + (f.*sin(theta_f)), Grd_my);

    # fishers move
    VR[VR.<.5] = 0.25 + (randn(P_C_n) * 0.1); # minimum speed
    CC_x = mod(CC[:,1] + (VR.*cos(theta_c)), Grd_mx);
    CC_y = mod(CC[:,2] + (VR.*sin(theta_c)), Grd_my);

    return CL_x,CL_y,CC_x,CC_y
end
