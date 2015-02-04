####!! Add processors
addprocs(2)

#### Add modules
@everywhere using Types
@everywhere import Iterators
@everywhere using PyCall
@everywhere using NPZ, Devectorize
@everywhere @pyimport rtree.index as pyrtree #R-tree python module
@everywhere unshift!(PyVector(pyimport("sys")["path"]), "") #require local folder 
@everywhere pypartition=pyimport(:random_partition) #Random partition for experiment

####!! Add functions and routines
require("Constants.jl");
require("sub_functions.jl");
require("sub_init.jl");
require("sub_routines.jl");
require("Experiments_core.jl");

####!! Experiement to parallelize
function do_nfisher_max()
    println("Optimal clique size vs all other sizes for 1 tauh and 1 taul")

    #! Parameters
	@everywhere function fnc_params()
		PRM.PC_lambda = 1.
		PRM.PS_n = int(round(.02/(pi*PRM.PF_sig^2/PRM.GRD_nx^2)));
		PRM.PC_n = 30; # number of fishers

		#! Get Tau_h and Tau_l from 2fisher case
		#! find these values by plotting fig3
		#! where there is MAXIMAL VOI
		data = npzread("./Data/Data_2fisher.npz")
		taul = data["taul"]
		tauh = data["tauh"]
		Hinf = data["H_inf"]
		Hnoinf = data["H_noinf"]
		Z = Hinf./Hnoinf
		S_p = 1 ./ taul
		C_q = PRM.PF_n ./ tauh

		id = find(Z.==maximum(Z))[1]
		PRM.PC_q = C_q[id];
		PRM.PS_p = S_p[id];

		PRM.PC_rp = fnc_optimal_PCrp();
		tauh = PRM.PF_n/PRM.PC_q
		taul = 1./PRM.PS_p
		taus = fnc_taus1()
		println("tau_h/tau_s: $(tauh/taus) tau_l/tau_s:  $(taul/taus)")

		#! Cliques
		occur = zeros(PRM.PC_n) #occurrence of a clique size
		Hrate = zeros(PRM.PC_n,2) #average catch rate per clique size
		Hdist = zeros(PRM.PC_n,2) #average CPUE per clique size
		TS = zeros(PRM.PC_n,2) #average search time per clique size
	end

	#! params
	fnc_params()
	
	#! run numerous clique trials
	@everywhere function fnc_trials(N)
		for i = 1:N # number of trials

			#random partition of PC_n
			part=pycall(pypartition["random_partition"],PyAny,PRM.PC_n)
			cliq=cell(length(part)) #composition of the cliques
			idx=1
			for p in 1:length(part)
				cliq[p]=idx:idx+part[p]-1
				idx+=part[p]
			end

			## Run
			println(i)
			school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
			fnc_init_sn(cons,FLAGS,cliq)
			FLAGS["ensemble"] = 1
			FLAGS["whilecond"] = 1
			make_ensemble(school,fish,cons,fishtree,EVENTS,FLAGS);

			### Running statistics (average Hrate within a clique)
			#for p in 1:length(part)
			#	#for each element in the partition
			#	occur[part[p]]+=1
			#	x=(cons.measure["Hrate"][cliq[p]])
			#	Hrate[part[p],:]+=[mean(x) mean(x)^2]
			#	x=(cons.measure["Hdist"][cliq[p]])
			#	Hdist[part[p],:]+=[mean(x) mean(x)^2]
			#	x=(cons.measure["Ts"][cliq[p]])
			#	TS[part[p],:]+=[mean(x) mean(x)^2]
			#end
		end
	end

	#! run it
	tic(); Y = pmap(fnc_trials,[1:2]); toc()

	tic(); fnc_trials(2); toc()

	#! Running statistics (average Hrate over cliques realizations)
    for p=1:PRM.PC_n
        if occur[p]>0
            Hrate[p,:]/=occur[p]
            Hdist[p,:]/=occur[p]
            TS[p,:]/=occur[p]
            #Variances
            Hrate[p,2]-=Hrate[p,1]^2
            Hdist[p,2]-=Hdist[p,1]^2
            TS[p,2]-=TS[p,1]^2
        end
    end

	#! SAVE
    npzwrite("./Data/Data_nfisher_max_p.npz", ["Hrate_mu"=>Hrate[:,1],"x"=>1,
            "Hdist_mu"=>Hdist[:,1],"TS_mu"=>TS[:,1],"Hrate_var"=>Hrate[:,2],
            "Hdist_var"=>Hdist[:,2],"TS_var"=>TS[:,2],"occur"=>occur])
end

