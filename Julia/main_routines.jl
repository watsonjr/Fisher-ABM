function sim_optimize()
ad = linspace(0,1,10);
TAU_MU = zeros(length(ad));
TAU_S2 = zeros(length(ad));
for i = 1:length(ad)

    ## modulate social network
    SN = ones(PC_n,PC_n) .* ad[i];
    for j = 1:PC_n; SN[j,j] = 1.; end;

    ## Initialize
    fish,cons,vars,tau = init_equilibrium();

    ## run model
    make_equilibrium(fish,cons,tau,vars,SN,0);

    ## calculate average encounter rate
    TAU_MU[i] = mean(tau.mu);
    TAU_S2[i] = mean(tau.s2);

end
return TAU_MU, TAU_S2
end








