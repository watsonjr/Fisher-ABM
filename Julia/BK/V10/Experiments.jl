
## simulate a simple scenario
function sim_simple()
sn = linspace(1e-6,1,10);   # types of prosociality
rp = 100; # number of repeats
CPUE = Array(Float64,length(sn),rp);
Tau  = Array(Float64,length(sn),rp);
for i = 1:length(sn)
	for j = 1:rp
		## modulate social network
		SN = ones(PC_n,PC_n) .* sn[i];
		for k = 1:PC_n; SN[k,k] = 1; end;

		## run model
		fish,cons,OUT = init_equilibrium();
		make_season(fish,cons,SN,0);

		## record
		CPUE[i,j] = mean(cons.cs ./ cons.Dist);
		Tau[i,j]  = mean(cons.Dist);
	end
	print(i/length(sn))
end
return CPUE, Tau
end

