
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
        DXY = [Dx[id,JJ][1] Dy[id,JJ][1]];
        DXY = DXY ./ sqrt(DXY[1].^2 + DXY[2].^2);

        if Dmin <= PC_c
            # probabilistic catch; if successful, then catch fish
            r =rand();if r < PC_q;KK = 1;else;KK = 0;end
        else
            KK = 0;
        end

    else # else I roam around randomly
        Dmin = 999; # flag, if you don't see anything
        JJ = 0; KK = 0;
        DXY = DXY + (randn(1,2).*PC_rn)
        DXY = DXY ./ sqrt(DXY[1].^2 + DXY[2].^2);
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
function fnc_move(CL,CC,DXY,D1)

    # clusters centers move randomly
    rn = rand(PCL_n);
    i = find(rn .<= PCL_p);
    CL_x = CL[:,1]; CL_y = CL[:,2];
    CL_x[i] = rand() * GRD_mx;
    CL_y[i] = rand() * GRD_my;

    # fishers move
    CC_x = mod(CC[:,1] + (DXY[:,1].*PC_v), GRD_mx);
    CC_y = mod(CC[:,2] + (DXY[:,2].*PC_v), GRD_my);

    return [CL_x CL_y],[CC_x CC_y], D1+1
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

## Running Stats
function fnc_stats(M,S,CPUE,TT)
     TT += 1;
     m  = M + ((CPUE-M)./TT);
     s  = S + ((CPUE-M).*(CPUE-m));
     return m,s,TT
end
