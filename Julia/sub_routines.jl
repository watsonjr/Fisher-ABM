


############## Run until encounter rates are stationary #############
function make_equilibrium(fish,cons,SN,ST)
TT = 1;
while sqrt(maximum(cons.ds2)) .> 0.0005
#for t = 1:2000

    ## Distances
    D,Dx,Dy = fnc_distance(fish.xy,cons.xy);

    ## Contact network from probabilistic social network
    fnc_contact(SN,cons.CN); # updates cons.CN

    ## Individual actions
    for i = 1:PC_n

        ## gather information
        #  Dmin: distance to nearest fish
        #  DDx,DDy: components of raw direction vector
        #  JJ: index of nearest fish
        (cons.Dmin[i],cons.DDx[i],cons.DDy[i],cons.JJ[i]) =
            fnc_information(D,Dx,Dy,cons.CN,i);

        ## direction
        #  ANG: direction vector
        #  VR:  travel speed (a function of Dmin)
        #  KK:  1/0; 1 if fish is caught, 0 if nothing caught
        (cons.ANG[i],cons.VR[i],cons.KK[i]) =
            fnc_direction(cons.Dmin[i],cons.DDx[i],cons.DDy[i],cons.ANG[i]);

    end

    ## Harvest
    (cons.H,fish.xy) = fnc_harvest(cons.KK,cons.JJ,cons.H,fish.xy);

    ## Move
    (fish.cl, cons.xy) = fnc_move(fish.cl,cons.xy,cons.ANG,cons.VR);

    ## Relocate fish (simulates movement)
    fish.xy = fnc_relocate(fish.cl,fish.xy);

    ## Fish growth
    fish.xy = fnc_growth(fish.xy,fish.cl);

    ## Statistics
    TT += 1; # ticker
    (cons.mu,cons.s2,cons.cs,cons.m,cons.s,cons.dmu,cons.ds2) =
        fnc_stats(cons.H,cons.cs,cons.s,cons.m,cons.mu,cons.s2,TT);

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





