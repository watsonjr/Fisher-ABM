


############## Run until encounter rates are stationary #############
function make_equilibrium(fish,cons,tau,vars,SN,ST)
while maximum(tau.dmu) > 0.01 || maximum(tau.ds2) > 1.0
    ## Distances
    D,Dx,Dy = fnc_distance(fish.xy,cons.xy);

    ## Information sharing
    for i = 1:PC_n

        ## gather information
        #  Dmin: distance to nearest fish
        #  DDx,DDy: components of raw direction vector
        #  JJ: index of nearest fish
        (vars.Dmin[i],vars.DDx[i],vars.DDy[i],vars.JJ[i])=fnc_information(D,Dx,Dy,SN,i);

        ## direction
        #  ANG: direction vector
        #  VR:  travel speed (a function of Dmin)
        #  KK:  1/0; 1 if fish is caught, 0 if nothing caught
        (vars.ANG[i],vars.VR[i],vars.KK[i]) =
            fnc_direction(vars.Dmin[i],vars.DDx[i],vars.DDy[i],vars.ANG[i]);

        ## harvest and update catch statistics (tau)
        (cons.H[i], fish.xy,
        tau.n[i],tau.s[i],tau.t[i],tau.mu[i],
        tau.M[i],tau.S[i],tau.s2[i],tau.dmu[i],tau.ds2[i]) =
        fnc_harvest_e(vars.KK[i],vars.JJ[i],cons.H[i],
                    fish.xy,fish.ci,fish.cl,
                    tau.n[i],tau.s[i],tau.t[i],tau.dmu[i],tau.mu[i],
                    tau.M[i],tau.S[i],tau.s2[i],tau.ds2[i]);

     end

    ## Move
    (fish.cl[:,1],fish.cl[:,2],cons.xy[:,1],cons.xy[:,2]) =
        fnc_move(fish.cl,cons.xy,vars.ANG,vars.VR);

    ## Relocate fish (simulates movement)
    fish.xy = fnc_relocate(fish.cl,fish.xy);

    ## Fish growth
    fish.xy = fnc_growth(fish.xy,fish.cl);
    

    ## Storage for plotting
    if ST == 1
        OUT.fish_xy = cat(3,OUT.fish_xy,fish.xy);
        OUT.cons_xy = cat(3,OUT.cons_xy,cons.xy);
    end
end
return
end





