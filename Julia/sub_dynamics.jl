

############## Run for a fixed length of time #############
function make_season_s(Fish_xy,Cons_xy,Cons_H,SN)
    for t = 1:P_Tend-1
        ## Distances
        D,Dx,Dy = fnc_distance(Fish_xy[:,:,t],Cons_xy[:,:,t]);

        ## Information sharing
        for i = 1:P_C_n

            ## gather information
            #  Dmin: distance to nearest fish
            #  DDx,DDy: components of raw direction vector
            #  JJ: index of nearest fish
            (Dmin[i], DDx[i], DDy[i], JJ[i]) = fnc_information(D,Dx,Dy,SN,i);

            (dmin1,ddx1,ddy1,jj1) = fnc_information(D,Dx,Dy,ones(P_C_n,P_C_n),i);
            (dmin2,ddx2,ddy2,jj2) = fnc_information(D,Dx,Dy,eye(P_C_n),i);

            ## direction
            #  ANG: direction vector
            #  VR:  travel speed (a function of Dmin)
            #  KK:  1/0; 1 if fish is caught, 0 if nothing caught
            (ANG[i], VR[i],KK[i]) =
                fnc_direction(Dmin[i],DDx[i],DDy[i],ANG[i]);

            ## harvest
            (Cons_H[i,t], Fish_xy[:,:,t]) =
                    fnc_harvest_s(KK[i],JJ[i],Cons_H[i,t],
                        Fish_xy[:,:,t],Fish_cl,Fclust_xy[:,:,t]);
         end

        ## Move
        (Fclust_xy[:,1,t+1],Fclust_xy[:,2,t+1],Cons_xy[:,1,t+1],Cons_xy[:,2,t+1]) =
            fnc_move(Fclust_xy[:,:,t],Cons_xy[:,:,t],ANG,VR);

        ## Relocate fish (simulates movement)
        Fish_xy[:,:,t+1] = fnc_relocate(Fclust_xy[:,:,t],Fish_xy[:,:,t]);
    end
    return Fish_xy,Cons_xy,Cons_H;
end



############## Run until encounter rates are stationary #############
function make_equilibrium(Fish_xy,Cons_xy,Cons_H,SN)
while max(Tau_dmu) > 0.01
    ## Distances
    D,Dx,Dy = fnc_distance(Fish_xy,Cons_xy);

    ## Information sharing
    for i = 1:P_C_n

        ## gather information
        #  Dmin: distance to nearest fish
        #  DDx,DDy: components of raw direction vector
        #  JJ: index of nearest fish
        (Dmin[i], DDx[i], DDy[i], JJ[i]) = fnc_information(D,Dx,Dy,SN,i);

        ## direction
        #  ANG: direction vector
        #  VR:  travel speed (a function of Dmin)
        #  KK:  1/0; 1 if fish is caught, 0 if nothing caught
        (ANG[i], VR[i],KK[i]) =
            fnc_direction(Dmin[i],DDx[i],DDy[i],ANG[i]);

        ## harvest
        (Cons_H[i], Fish_xy,
        Tau_n[i],Tau_s[i],Tau_t[i],Tau_dmu[i],Tau_mu[i]) =
                fnc_harvest_e(KK[i],JJ[i],Cons_H[i],
                    Fish_xy,Fish_cl,Fclust_xy,
                    Tau_n[i],Tau_s[i],Tau_t[i],Tau_dmu[i],Tau_mu[i]);

     end

    ## Move
    (Fclust_xy[:,1],Fclust_xy[:,2],Cons_xy[:,1],Cons_xy[:,2]) =
        fnc_move(Fclust_xy,Cons_xy,ANG,VR);

    ## Relocate fish (simulates movement)
    Fish_xy = fnc_relocate(Fclust_xy,Fish_xy);
end
return Tau_mu
end



