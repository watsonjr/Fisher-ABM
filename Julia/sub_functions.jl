
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


#### FISH RELOCATION
function fnc_relocate(Cl,FF)
    i = int(ceil(rand(1)*(PF_m.*PF_n))); # choose a number of fish to relocate
    j = shuffle([1:PF_n]); # shuffle for randomness
    j = j[1:i[1]]; # choose the select few
    k = ceil(rand(i[1],1) * PCL_n); # cluster center to move to
    a = repmat([GRD_mx GRD_my],length(j),1); # for periodic boundaries
    F = mod(Cl[k[:],:] + randn(i[1],2).*PF_dx, a); # move them
    FF[j,:] = F; # add back
    return FF
end


#### DISTANCES
function fnc_distance(FF,CC)
    ## x,y, distances
    Dx = FF[:,1]' .- CC[:,1];
    Dy = FF[:,2]' .- CC[:,2];

    ## periodic boundaries
    i = (abs(Dx).>(GRD_mx/2)) + (abs(Dy).>(GRD_my/2));
    Dx[i.>1] = -sign(Dx[i.>1]).*(GRD_mx - abs(Dx[i.>1]));
    Dy[i.>1] = -sign(Dy[i.>1]).*(GRD_my - abs(Dy[i.>1]));

    ## Distance
    D = sqrt(Dx.^2 + Dy.^2);
    return D,Dx,Dy
end


#### INFORMATION
function fnc_information(D,Dx,Dy,CN,id)
    # Who you will get info from
    j = find(CN[id,:].==1);

    # get information from friends and yourself
    DDi = 0;
    for k = 1:length(j)
        d = D[j[k],:];
        kk = find(d.<PC_ff); # find only those in fish finder view
        DDi = hcat(DDi,kk');
    end
    DDi = DDi[2:end];
    DD  = D[id,DDi];

    # find nearest fish
    if isempty(DD)==0 # if I see something
        Dmin = minimum(DD);
        II = find(DD .== Dmin);
        JJ = DDi[II][1];
        DDx = Dx[id,JJ][1];
        DDy = Dy[id,JJ][1];

    else # else I roam around randomly
        Dmin = 999; # flag, if you don't see anything
        JJ = 0;
        DDx  = 0; rand();
        DDy  = 0; rand();
    end

    return Dmin,DDx,DDy,JJ
end


#### DIRECTION and SPED
function fnc_direction(Dmin,DDx,DDy,ANG_in)
    if (Dmin <= PC_c && Dmin <= PC_ff) # if very near - go fishing
        DXY = [DDx DDy] ./ Dmin;
        ANG_out  = DXY[:,2] ./ DXY[:,1];
        vr = 0.001; # no speed, stay where you are

        # probabilistic catch
        r =rand();
        if r < PC_pr # if successful, then catch fish
            KK  = 1;
        else
            KK = 0;
        end

    elseif (Dmin > PC_c && Dmin <= PC_ff) # if nearish move towards slowly
        DXY  = [DDx DDy] ./ Dmin;
        ANG_out  = DXY[:,2] ./ DXY[:,1];
        vr   = (Dmin ./ PC_ff) .* PC_v;
        KK   = 0;

    elseif (Dmin > PC_c && Dmin > PC_ff && Dmin < 999) # if not nearish move towards fast
        DXY  = [DDx DDy] ./ Dmin;
        ANG_out  = DXY[:,2] ./ DXY[:,1];
        vr   = PC_v;
        KK   = 0;

    else  # if not so near, random walk
        ANG_out = ANG_in + (PC_ang * rand() - (PC_ang/2));
        vr   = PC_v;
        KK   = 0;

    end
    return ANG_out[1],vr,KK
end


#### HARVEST for a season
function fnc_harvest(KK,JJ,CC,FF);
    II    = KK.*JJ
    IIu   = unique(II);

    for i in IIu
        j = find(II.==i);
        CC[j] =  (1 / length(j)) * KK[j];
    end
    FF[IIu[IIu.>0],:] = NaN;
    return CC,FF
end

#function fnc_harvest(KK,JJ,CC,Tau_n,Tau_t,
#                       Tau_s,Tau_mu,Tau_S,Tau_M,Tau_s2,Tau_dmu,Tau_ds2,FF)
#    if II != 0 # if the fish is caught
#
#        # record that is it
#        CC = 1; # catch
#
#        # update waiting time mean
#        tau_n   = Tau_n + 1; # increment catch event counter
#        tau_s   = Tau_s + Tau_t; # accumulate waiting times
#        tau_mu  = tau_s / tau_n;
#        tau_dmu = abs(Tau_mu - tau_mu);
#        tau_t   = 1; # reset wait time counter
#
#        # update waiting time variance
#        tau_M   = Tau_M + ((Tau_t - Tau_M) / tau_n);
#        tau_S   = Tau_S + ((Tau_t-Tau_M) * (Tau_t-tau_M));
#        tau_s2  = tau_S / (tau_n-1);
#        tau_ds2 = abs(Tau_s2 - tau_s2);
#
#        # kill fish
#        FF[II,:] = NaN;
#
#    else
#
#        CC = 0; # no catch
#        tau_t  = Tau_t + 1; # increment wait time counter
#
#        # update waiting time mean
#        tau_n  = Tau_n;
#        tau_s  = Tau_s;
#        tau_mu = Tau_mu;
#
#        # update waiting time variance
#        tau_S  = Tau_S;
#        tau_M  = Tau_M;
#        tau_s2 = Tau_s2;
#        tau_dmu= Tau_dmu;
#        tau_ds2= Tau_ds2;
#    end
#    return CC,tau_n,tau_t,tau_s,tau_mu,tau_S,tau_M,tau_s2,tau_dmu,tau_ds2,FF
#end


#### MOVE
function fnc_move(CL,CC,theta_c,VR)

    # clusters centers move randomly
    rn = rand(PCL_n);
    i = find(rn .<= PCL_p);
    CL_x = CL[:,1]; CL_y = CL[:,2];
    CL_x[i] = rand() * GRD_mx;
    CL_y[i] = rand() * GRD_my;

    # fishers move
    CC_x = mod(CC[:,1] + (VR.*cos(theta_c)), GRD_mx);
    CC_y = mod(CC[:,2] + (VR.*sin(theta_c)), GRD_my);

    return [CL_x CL_y],[CC_x CC_y]
end


#### Fish Growth
function fnc_growth(FF,CL);
    i = find(isnan(FF[:,1])); # find nans (i.e. dead fish)

    if isempty(i)==false
        rn = rand()
        j = shuffle(i);
        if rn <= PF_g
            k = ceil(rand() * PCL_n); # cluster center to move to
            F = mod(CL[k,:] + randn().*PF_dx, [GRD_mx GRD_my]); # move them
            FF[j[1],:] = F;
        end
    end
return FF
end
