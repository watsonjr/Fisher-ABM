
#### Run over a whole season
function make_season(Fish_xy,Cons_xy,Cons_H)

    for t = 1:P_Tend-1

        ## Distances
        D,Dx,Dy = fnc_distance(Fish_xy[:,:,t],Cons_xy[:,:,t]);

        ## Information sharing
        for i = 1:P_C_n

            ## gather information
            #  Dmin: distance to nearest fish
            #  DDx,DDy: components of raw direction vector
            #  JJ: index of nearest fish
            (Dmin[i], DDx[i], DDy[i], JJ[i]) = fnc_information(D,Dx,Dy,i);

            ## direction
            #  ANG: direction vector
            #  VR:  travel speed (a function of Dmin)
            #  KK:  1/0; 1 if fish is caught, 0 if nothing caught
            (ANG[i], VR[i],KK[i]) =
                fnc_direction(Dmin[i],DDx[i],DDy[i],ANG[i]);

            ## harvest
            (Cons_H[i,t], Fish_xy[:,:,t]) =
                fnc_harvest(KK[i],JJ[i],Cons_H[i,t],
                            Fish_xy[:,:,t],Fish_cl,Fclust_xy[:,:,t]);

         end

        ## Move
        (Fclust_xy[:,1,t+1],Fclust_xy[:,2,t+1],Cons_xy[:,1,t+1],Cons_xy[:,2,t+1]) =
            fnc_move(Fclust_xy[:,:,t],Cons_xy[:,:,t],ANG,VR);

        ## Relocate fish (simulates movement)
        Fish_xy[:,:,t+1] = fnc_relocate(Fclust_xy[:,:,t],Fish_xy[:,:,t]);
    end

    return Fish_xy,Cons_xy,Cons_H
end
