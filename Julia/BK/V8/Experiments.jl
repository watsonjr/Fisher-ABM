
##### RUN THE OPTIMIZATION MODEL
function sim_optimize()
ad = linspace(1e-6,1,10);   # types of prosociality
pr = linspace(0.01,1,10);
MU = zeros(length(pr),length(ad));
S2 = zeros(length(pr),length(ad));
global PC_pr
for i = 1:length(ad)
    for j = 1:length(pr)
        ## modulate fish distribution
        PC_pr = pr[j];

        ## modulate social network
        SN = ones(PC_n,PC_n) .* ad[i];
        for k = 1:PC_n; SN[k,k] = 1; end;

        ## run for equilibriumm
        M,S,TT = make_equilibrium(SN);

        ## Store
        MU[i,j] = mean(M);
        S2[i,j] = mean(S ./ (TT-1));
    end
    println(i/length(ad)," | ",j/length(pr));
end
return MU,S2
end



##### RUN THE EVOLUTIONARY MODEL
function sim_evolve()
# setup
P_tend = 200;
M1 = zeros(PC_n);
S1 = zeros(PC_n);
SN = ones(PC_n,PC_n) .* 0.25;
for j = 1:PC_n; SN[j,j] = 1.; end;
OUT_SN      = zeros(PC_n,PC_n,P_tend);
OUT_CPUE_mu = zeros(PC_n,P_tend);
OUT_CPUE_s2 = zeros(PC_n,P_tend);

# run
for t = 1:P_tend # need to replace this with a while loop

    ## Initialize
    fish,cons = init_equilibrium();
    M2,S2 = make_equilibrium(SN);

    ## calculate change in fitness
    DF = M2 - M1;

    ## update social preference
    cons.S = fnc_makebreak(DF,cons.S);

    ## update social network
    SN = fnc_update_SN(SN,cons.S);

    ## update running stats
    M1 = M2;
    S1 = S2;

    ## Save
    OUT_SN[:,:,t] = SN;
    OUT_CPUE_mu[:,t] = M2;
    OUT_CPUE_s2[:,t] = S2;
    t = t+1;
    println(t);

end
npzwrite("./Data/Data_evo_SN.npy", OUT_SN)
npzwrite("./Data/Data_evo_cpue_mu.npy", OUT_CPUE_mu)
npzwrite("./Data/Data_evo_cpue_s2.npy", OUT_CPUE_s2)


return OUT_SN, OUT_CPUE_mu, OUT_CPUE_s2
end







