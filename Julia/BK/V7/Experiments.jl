
##### RUN THE OPTIMIZATION MODEL
function sim_optimize()

ad = linspace(1e-6,1,10);
MU = zeros(PC_n,length(ad));
S2 = zeros(PC_n,length(ad));
for i = 1:length(ad)

    ## modulate social network
    SN = ones(PC_n,PC_n) .* ad[i];
    for j = 1:PC_n; SN[j,j] = 1.; end;

    ## Loop multiple seasons
    mu = zeros(PC_n,10);
    s2 = zeros(PC_n,10);
    for j = 1:10

        ## Initialize
        fish,cons = init_equilibrium();

        ## run model
        make_equilibrium(fish,cons,SN,0);

        ## save
        mu[:,j] = cons.cs;

    end

    ## calculate CPUE stats
    MU[:,i] = mean(mu,2);
    S2[:,i] = var(mu,2);

    ## ticker
    println(i/length(ad));
end
S = ad;
npzwrite("./Data/Data_opt_cpue_mu.npy", MU)
npzwrite("./Data/Data_opt_cpue_s2.npy", S2)
npzwrite("./Data/Data_opt_S.npy", S)


return MU, S2, ad
end


##### RUN THE EVOLUTIONARY MODEL
function sim_evolve()
# setup
P_tend = 600;
MU_1 = zeros(PC_n);
MU_2 = zeros(PC_n);
SN   = ones(PC_n,PC_n) .* 0.25;
for j = 1:PC_n; SN[j,j] = 1.; end;
OUT_SN      = zeros(PC_n,PC_n,P_tend);
OUT_CPUE_mu = zeros(PC_n,P_tend);
OUT_CPUE_s2 = zeros(PC_n,P_tend);

# run
for t = 1:P_tend # need to replace this with a while loop

    ## Initialize
    fish,cons = init_equilibrium();

    ## Record current fitness
    MU_1 = cons.mu;

    ## run model
    make_equilibrium(fish,cons,SN,0);

    ## calculate average encounter rate
    MU_2 = cons.mu;

    ## calculate change in fitness
    DF = MU_2 - MU_1;

    ## update social preference
    cons.S = fnc_makebreak(DF,cons.S);

    ## update social network
    SN = fnc_update_SN(SN,cons.S);

    ## Save
    OUT_SN[:,:,t] = SN;
    OUT_CPUE_mu[:,t] = cons.mu;
    OUT_CPUE_s2[:,t] = cons.s2;
    t = t+1;
    println(t);

end
return OUT_SN, OUT_CPUE_mu, OUT_CPUE_s2
end







