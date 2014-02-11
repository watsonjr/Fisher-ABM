

#### Preamble
include("make_parameters.jl");
using NPZ
using Distance

#### Fish movement
for t = 1:P_Tend-1

    ######## FISH
    ###### Grid fish density
    #ii = int(round(Fish[:,:,1]));
    #jj,mm,kk = hist2d(ii,Grd_x[1:2:end],Grd_y[1:2:end]);

    ###### Fish demographics
    i = int(ceil(rand(1)*(F_n./100))); # choose a number of fish to kill
    j = shuffle([1:F_n]);
    j = j[1:i[1]];

    ## Relocate fish
    k = ceil(rand(i[1],1) * Fc_n); # cl center to move to
    a = repmat([Grd_mx Grd_my],length(j),1);
    Fish[j,:,t] = mod(Fclust[k[:],:,t] + randn(i[1],2).*F_dx, a);

    ######## FISHERS
    #### Find nearest fish
    ## x,y, distances
    Dx = Fish[:,1,t]' .- Cons[:,1,t];
    Dy = Fish[:,2,t]' .- Cons[:,2,t];

    ## periodic boundaries
    i = (abs(Dx).>(Grd_mx/2)) + (abs(Dy).>(Grd_my/2));
    Dx[i.>1] = -sign(Dx[i.>1]).*(Grd_mx - abs(Dx[i.>1]));
    Dy[i.>1] = -sign(Dy[i.>1]).*(Grd_my - abs(Dy[i.>1]));

    ## Distance
    D = sqrt(Dx.^2 + Dy.^2);

    ## Find fish in fishfinder
    for i = 1:C_n
        d = min(D[i,:]);
        j = find(D[i,:].==d);

        if d <= P_ddx && d >= P_catch # if nearish move towards
            DXY[i,:] = [Dx[i,j] Dy[i,j]] ./ D[i,j];
            Cnr[i] = D[i,j][1] ./ P_ddx;
        elseif d < P_catch # if very near - go fishing
            r =rand();
            if r < P_prob # if successful
                Harv[i,t] = 1;
                k = ceil(rand() * Fc_n); # cl center to move to
                Fish[j,:,t] = mod(Fclust[k,:,t]+randn().*F_dx,[Grd_mx Grd_my]);
            end
        else # if not so near random walk
            dx = rand();
            dy = rand();
            dd = sqrt(dx.^2 + dy.^2);
            DXY[i,:] = [dx dy] ./ dd;
            Cnr[i] = 1;
        end
    end

    #for i = 1:C_n
    #    cx = Cons[i,1,t]; cy = Cons[i,2,t];
    #    fx = Fish[:,1,t]; fy = Fish[:,2,t];
    #    sx = sign(cx-fx); sy = sign(cy-fy);
    #    dx = mod(cx-fx,Grd_mx/2).*sx; dy = mod(cy-fy,Grd_my/2).*sy;
    #    j  = find(((abs(dx).<P_ddx) + (abs(dy).<P_ddx)).==2);#index fish in fishfinder

    #    if isempty(j) == false
    #        D  = sqrt((cx - Fish[j,1,t]).^2 + (cy-Fish[j,2,t]).^2);
    #        k  = find(D.==min(D)); # index of nearest fish
    #        DXY[i,:] = [[dx[j[k]] dy[j[k]]] ./ D[k]]; # norm direction to nearest fish
    #        Cnr[i,1] = D[k][1] ./ P_ddx; # distance to nearest 'visible' fish
    #    else
    #        dxy = [rand() rand()];
    #        DXY[i,:] = dxy ./ sqrt(sum(dxy.^2));
    #        Cnr[i,1] = 1; # flag, no fish near
    #    end
    #end

    ######## MOVE
    ###### Move clusters centers and fish
    theta_f = rand(Fc_n,1)*2*pi;
    f = rand(Fc_n,1).^(-1/Fc_alp);
    Fclust[:,1,t+1] = mod(Fclust[:,1,t] + (f.*cos(theta_f)), Grd_mx);
    Fclust[:,2,t+1] = mod(Fclust[:,2,t] + (f.*sin(theta_f)), Grd_my);

    ###### Move fishers
    Cons[:,1,t+1] = mod(Cons[:,1,t] + (Cnr.*C_v.*DXY[:,1]), Grd_mx);
    Cons[:,2,t+1] = mod(Cons[:,2,t] + (Cnr.*C_v.*DXY[:,2]), Grd_my);

    ######## update
    Fish[:,:,t+1] = Fish[:,:,t];
end


#### Save
npzwrite("./Data/Data_fish.npy", Fish)
npzwrite("./Data/Data_fishers.npy", Cons)
npzwrite("./Data/Data_fclust.npy", Fclust)


