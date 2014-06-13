
#### Run a season
function make_season(fish,cons,SN,ST)
while minimum(cons.cs) .< PF_n

    ## Distances
    D,Dx,Dy = fnc_distance(fish.xy,cons.xy);

    ## Contact network from probabilistic social network
    fnc_contact(SN,cons.CN); # updates cons.CN

    ## Individual actions
    for i = 1:PC_n

        ## gather information
        (cons.Dmin[i],cons.DXY[i,:],cons.JJ[i],cons.KK[i]) =
            fnc_information(D,Dx,Dy,cons.DXY[i,:],cons.CN,i);

    end

    ## Harvest
    (cons.H,cons.cs,fish.xy) = fnc_harvest(cons.KK,cons.JJ,cons.H,cons.cs,fish.xy);

    ## Move
    (fish.cl,cons.xy,cons.Dist) = fnc_move(fish.cl,cons.xy,cons.DXY,cons.Dist);

    ## Relocate fish (simulates movement)
    fish.xy = fnc_relocate(fish.cl,fish.xy);

    ## Fish growth
    fish.xy = fnc_growth(fish.xy,fish.cl);

    ## Storage for plotting
    if ST == 1
        OUT.fish_xy = cat(3,OUT.fish_xy,fish.xy);
        OUT.cons_xy = cat(3,OUT.cons_xy,cons.xy);
        OUT.clus_xy = cat(3,OUT.clus_xy,fish.cl);
        OUT.cons_H = cat(2,OUT.cons_H,cons.H);
    end
end
return
end


#### Run a season
function make_equilibrium(SN)
    ## setup running stats
    TT = 1;
    M  = zeros(PC_n); # running mean CPUE
    S  = zeros(PC_n); # running variance CPUE

    ## Run initial season
    fish,cons = init_equilibrium();
    make_season(fish,cons,SN,0);
    M = cons.cs ./ cons.Dist; # initial CPUe

    ## Run second season
    fish,cons = init_equilibrium();
    make_season(fish,cons,SN,0);
    M,S,TT = fnc_stats(M,S,cons.cs./cons.Dist,TT);
    df = M;

    ## Loop multiple seasons
    while maximum(df./M) .> 0.1
        ## run model
        fish,cons = init_equilibrium();
        make_season(fish,cons,SN,0);

        ## calculate running stats
        m,S,TT = fnc_stats(M,S,cons.cs./cons.Dist,TT);

        ## update running stats
        df = abs(M-m); M = m;
    end
    return M,S,TT
end


