
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
function fnc_information(D,Dx,Dy,DXY,CN,id)
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
        DDx  = DXY[1]; rand();
        DDy  = DXY[2]; rand();
    end

    return Dmin,[DDx DDy],JJ
end


#### DIRECTION and SPED
function fnc_direction(Dmin,DXY)
    if (Dmin <= PC_c && Dmin <= PC_ff && Dmin < 999) # if very near - go fishing
        DXY = DXY ./ Dmin;
        vr = 0.001; # no speed, stay where you are

        # probabilistic catch; if successful, then catch fish
        r =rand();if r < PC_pr;KK = 1;else;KK = 0;end

    elseif (Dmin > PC_c && Dmin <= PC_ff && Dmin < 999) # if nearish move towards slowly
        DXY  = DXY ./ Dmin;
        vr   = (Dmin ./ PC_ff) .* PC_v;
        vr   = vr + (vr * rand() * 0.1);
        KK   = 0;

    elseif (Dmin > PC_c && Dmin > PC_ff && Dmin < 999) # if not nearish move towards fast
        DXY  = DXY ./ Dmin;
        vr   = PC_v + (PC_v * rand() * 0.1);
        KK   = 0;

    else  # if not so near, random walk
        DXY = DXY + (randn(1,2).*PC_rn)
        DXY = DXY ./ sqrt(DXY[1].^2 + DXY[2].^2);
        vr   = PC_v + (PC_v * rand() * 0.1);
        KK   = 0;

    end
    return DXY,vr,KK
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


#### STATISTICS
function fnc_stats(H,CS,m,s,Dist,TT);
#H = cons.H; # harvest this time step
#CS = cons.cs; # cumulative harvest
#m = cons.m; # component of variance calculation
#s = cons.s; # component of variance caluclation
#Dist = cons.Dist; # current distance traveled
#TT = TT; # current time step

    cs   = CS + H; # cumulative catch
    cpue = cs ./ Dist; # cumulative cpue

    M    = m + ((cpue-m)./TT); # running mean cpue
    S    = s + ((cpue-m).*(cpue-M)); # running variance in cpue

    s2   = s ./ (TT-1);
    S2   = S ./ (TT-1);

    dmu = abs(m-M);
    ds2 = abs(s2-S2);

    return cs,M,S,S2,dmu,ds2
end


#### MOVE
function fnc_move(CL,CC,DXY,D1,VR)

    # clusters centers move randomly
    rn = rand(PCL_n);
    i = find(rn .<= PCL_p);
    CL_x = CL[:,1]; CL_y = CL[:,2];
    CL_x[i] = rand() * GRD_mx;
    CL_y[i] = rand() * GRD_my;

    # fishers move
    #dx = (VR.*cos(theta_c));
    #dy = (VR.*sin(theta_c));
    D2  = sqrt(DXY[:,1].^2 + DXY[:,2].^2);
    CC_x = mod(CC[:,1] + (DXY[:,1].*VR), GRD_mx);
    CC_y = mod(CC[:,2] + (DXY[:,2].*VR), GRD_my);

    return [CL_x CL_y],[CC_x CC_y], D1 + D2
end


#### Fish Growth
function fnc_growth(FF,CL);
    i = find(isnan(FF[:,1])); # find nans (i.e. dead fish)

    if isempty(i)==false
        rn = rand()
        j = shuffle(i);
        if rn <= PF_g
            k = ceil(rand() * PCL_n); # cluster center to move to
            Fx = mod(CL[k,1] + randn().*PF_dx, GRD_mx); # move them
            Fy = mod(CL[k,2] + randn().*PF_dx, GRD_my); # move them
            FF[j[1],:] = [Fx Fy];
        end
    end
return FF
end


#### Update social preference (make/break friendships);
function fnc_makebreak(DF,S)
    for i = 1:PC_n
        if DF[i] .> 0;
            # stick
            S[i] = S[i] .* 1;
        elseif DF[i] .< 0
            # twist
            S[i] = S[i] .* -1;
        end
    end
    return S
end

#### Update social network
function fnc_update_SN(SN,S);
    for i = 1:PC_n

        # select a fisher to update make/break
        sn = SN[i,:];

        if S[i] == 1
            sn[i] = 0; # you don't change your relationship with yourself
            sn[sn.==1] = 0; # cant improve bbfs
            sn  = sn ./ sum(sn);
            csn = cumsum(sn,2);
            rn  = rand();
            Del = csn-rn;
        else
            sn[i] = 0; # don't hate yourself
            sn  = 1 ./ sn; # inverse
            sn[isinf(sn).==1] = 0; # don't hate your worst enemies more
            sn  = sn ./ sum(sn);
            csn = cumsum(sn,2);
            rn  = rand();
            Del = csn-rn;
        end

        # pick a fisher
        j = find(Del.>0); # pick one

        # change social network
        if isempty(j) == 0
            j = j[1];
            SN[i,j] = SN[i,j] + (S[i] * (0.05*rand()));
        end
    end

    # make sure they're a probability
    SN[SN.<=0] = 1e-6;
    SN[SN.>1]  = 1.0;

    return SN
end

