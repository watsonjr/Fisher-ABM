

#### Preamble
include("make_parameters.jl");
using NPZ
using Distance

#### Fish_xy movement
for t = 1:P_Tend-1

    ######## FISH
    ## Relocate fish (simulates movement)
    i = int(ceil(rand(1)*(P_F_m.*P_F_n))); # choose a number of fish to kill
    j = shuffle([1:P_F_n]);
    j = j[1:i[1]];
    k = ceil(rand(i[1],1) * P_Cl_n); # cl center to move to
    a = repmat([Grd_mx Grd_my],length(j),1);
    Fish_xy[j,:,t] = mod(Fclust_xy[k[:],:,t] + randn(i[1],2).*P_F_dx, a);

    ######## FISHERS
    #### Find nearest fish
    ## x,y, distances
    Dx = Fish_xy[:,1,t]' .- Cons_xy[:,1,t];
    Dy = Fish_xy[:,2,t]' .- Cons_xy[:,2,t];

    ## periodic boundaries
    i = (abs(Dx).>(Grd_mx/2)) + (abs(Dy).>(Grd_my/2));
    Dx[i.>1] = -sign(Dx[i.>1]).*(Grd_mx - abs(Dx[i.>1]));
    Dy[i.>1] = -sign(Dy[i.>1]).*(Grd_my - abs(Dy[i.>1]));

    ## Distance
    D = sqrt(Dx.^2 + Dy.^2);

    ## Information sharing
    for i = 1:P_C_n
        j = find(P_C_Sn[i,:] .== 1); # who is in my social network

        # get information from friends and yourself
        DD = 0; DDx = 0; DDy = 0;
        for k = 1:length(j)
            d = D[j[k],:];
            DD  = hcat(DD,D[j[k],:]); ### Choose only those fish in finder
            DDx = hcat(DDx,Dx[j[k],:]);
            DDy = hcat(DDy,Dy[j[k],:]);
        end
        DD  = DD[2:end];
        DDx = DDx[2:end];
        DDy = DDy[2:end];

        # find nearest fish ### Could change this to greatest density of fish
        dmin = min(DD);
        j = find(DD .== dmin);

        # get direction vector
        if dmin <= P_ddx && dmin >= P_catch # if nearish move towards
            P_DXY[i,:] = [Dx[i,j] Dy[i,j]] ./ D[i,j];
            P_vr[i] = D[i,j][1] ./ P_ddx;
        elseif dmin < P_catch # if very near - go fishing
            r =rand();
            if r < P_prob # if successful
                Cons_H[i,t] = 1;
                k = ceil(rand() * P_Cl_n); # cl center to move to
                Fish_xy[j,:,t] = mod(Fclust_xy[k,:,t]+randn().*P_F_dx,[Grd_mx Grd_my]);
            end
        else # if not so near random walk
            dx = rand();
            dy = rand();
            dd = sqrt(dx.^2 + dy.^2);
            P_DXY[i,:] = [dx dy] ./ dd;
            P_vr[i] = 1;
        end
    end

    ######## MOVE
    ###### Move clusters centers and fish
    theta_f = rand(P_Cl_n,1)*2*pi;
    f = rand(P_Cl_n,1).^(-1/P_Cl_al);
    Fclust_xy[:,1,t+1] = mod(Fclust_xy[:,1,t] + (f.*cos(theta_f)), Grd_mx);
    Fclust_xy[:,2,t+1] = mod(Fclust_xy[:,2,t] + (f.*sin(theta_f)), Grd_my);

    ###### Move fishers
    Cons_xy[:,1,t+1] = mod(Cons_xy[:,1,t] + (P_vr.*P_C_v.*P_DXY[:,1]), Grd_mx);
    Cons_xy[:,2,t+1] = mod(Cons_xy[:,2,t] + (P_vr.*P_C_v.*P_DXY[:,2]), Grd_my);

    ######## update
    Fish_xy[:,:,t+1] = Fish_xy[:,:,t];
end


#### Save
npzwrite("./Data/Data_fish.npy", Fish_xy)
npzwrite("./Data/Data_fishers.npy", Cons_xy)
npzwrite("./Data/Data_fclust.npy", Fclust_xy)

#### Useful Julia code
###### Grid fish density
#ii = int(round(Fish_xy[:,:,1]));
#jj,mm,kk = hist2d(ii,Grd_x[1:2:end],Grd_y[1:2:end]);

