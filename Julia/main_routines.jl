
##### RUN THE OPTIMIZATION MODEL
function sim_optimize()
ad = linspace(0,1,10);
TAU_MU = zeros(length(ad));
TAU_S2 = zeros(length(ad));
for i = 1:length(ad)

    ## modulate social network
    SN = ones(PC_n,PC_n) .* ad[i];
    for j = 1:PC_n; SN[j,j] = 1.; end;

    ## Initialize
    fish,cons = init_equilibrium();

    ## run model
    make_equilibrium(fish,cons,SN,0);

    ## calculate average encounter rate
    TAU_MU[i] = mean(cons.mu);
    TAU_S2[i] = mean(cons.s2);

end
return TAU_MU, TAU_S2
end


##### RUN THE EVOLUTIONARY MODEL
function sim_evolve()
MU_1 = zeros(PC_n);
MU_2 = zeros(PC_n);
P_tend = 1000;
OUT_SN = zeros(size(SN,1),size(SN,2),P_tend);
for t = 1:P_tend # need to replace this with a while loop

    ## Record current fitness
    MU_1 = cons.mu;

    ## Initialize
    fish,cons = init_equilibrium();

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
    t = t+1;
    println(t);

end
return OUT_SN
end







