
##########################!!!!! SIMULATION EXPERIMENTS !!!!!!!!################
##################### BASIC TESTS
#! Timing test: ... this does what
function do_timingtest(fsave=true)
    println("Timing test")

	## Choose parameters
    PRM.PC_lambda = 1. # strength of social network tie
    PRM.PC_n = 2; # number of fishers

    ## Initialize
    school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
    fnc_init_sn(cons,FLAGS) 

    ## Force SN
    cons.SN = ones(PRM.PC_n,PRM.PC_n);
    #cons.SN = eye(PRM.PC_n,PRM.PC_n);

    ## Choose flags
    FLAGS["save"] = fsave
    FLAGS["whilecond"] = 1 # short run
    @time  make_ensemble(school,fish,cons,fishtree,EVENTS,FLAGS,OUT);
end


##################### 1 SCHOOl 1 FISHER
function do_1fisher()
    ######  Expected search time (T_s^*)
    println("Expected search time")
    NN = 20; #! number of replicates
    RP = linspace(0.1,.9,20); #! probability of turning per unit time
    result = cell(length(RP),NN); #! array for results collection

    ##! Parameters
    PRM.PC_n = 1

    #! Setup save
    S_rp = zeros(NN); S_Ts = zeros(length(RP),NN) ; S_Ts1 = zeros(NN);

    for i = 1:length(RP) # Iterate over different turning probs
    	for j = 1:NN # run sim a number of times 
    		## initialize
			print(i/length(RP)," | ", j/NN,"\n")
			PRM.PC_rp = RP[i];
			school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
			fnc_init_sn(cons,FLAGS)

			## flags
			FLAGS["ensemble"] = 1; # only one season
			FLAGS["whilecond"] = 1; # stop when harvest rate has converged
			
			## Run
			make_ensemble(school,fish,cons,fishtree,EVENTS,FLAGS);

			# save results
			#result[i,j] = RP[i],mean(cons.measure["Ts"]),
			#			mean(cons.measure["Tv"]);
			@inbounds S_Ts[i,j] = mean(cons.measure["Ts"])
		end
		S_rp[i] = RP[i];
		S_Ts1[i] = fnc_taus1();
    end
    #save_results(result,R"PC_rp \tau_s^R std","Data_1fisher_1school_a",PRM)
    npzwrite("./Data/Data_1fisher.npz", ["PC_rp"=>S_rp,"x"=>1,"Ts"=>S_Ts,"Ts1"=>S_Ts1])

    #! check
	#PRM.PC_rp = fnc_optimal_PCrp()
	#taus_min = fnc_taus1()
end


#function do_1fisher_1school_b()
#    ###### T_s^* as a function of C_f and F_sig
#    SIG = linspace(2,GRD_mx/10,25); # F_sigma (radius of school)
#    FF  = linspace(2,GRD_mx/10,25); # fish finder radius
#    result = cell(size(SIG,1),size(FF,1));
#
#    ##! Parameters
#    PRM.PC_n=1
#
#    ## Run
#    for i = 1:length(SIG)
#        for j = 1:length(FF)
#            print(i," ",j,"\n")
#
#            ## Initialize
#            PRM.PF_sig = SIG[i];
#            PRM.PC_f   = FF[j];
#            PRM.PC_rp  = fnc_optimal_PCrp();
#            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
#            fnc_init_sn(cons,FLAGS)
#
#            ## Flags
#			FLAGS["ensemble"] = 5
#			FLAGS["whilecond"] = 1
#
#			## Run and save
#            make_ensemble(school,fish,cons,fishtree,EVENTS,FLAGS);
#            result[i,j] =SIG[i],FF[j], mean(cons.measure["Ts"]), 
#            	mean(cons.measure["Tv"]),PRM.PC_rp;
#        end
#    end
#    save_results(result,R"PF_sig PC_f \tau_s^R std PC_rp",
#    	"Data_1fisher_1school_b",PRM)
#end


##################### 1 SCHOOL 2 FISHERS
##! Test effect of friendship, on Tau_s and CPUE, 
##! as Tau_l and Tau_h are varied
function do_2fisher()
    ###### Test performance as a function of tau_l and tau_h
    println("Performance as a function of tau_l and tau_h")
    println("Number of schools: $(PRM.PS_n)")

    ## basic values of C_f and F_sig
    PRM.PC_n = 2; #2 fishers!
    #PRM.PS_n = int(round(.02/(pi*PRM.PF_sig^2/PRM.GRD_nx^2))); #fish area
    
    #! Choose resolution of Tau_l and Tau_h
    JUMP   = logspace(log10(0.001),log10(.1),25); ## PS_p
    CATCH  = logspace(log10(0.01),log10(2),25); ## PC_q
    result_inf   = cell(size(JUMP,1),size(CATCH,1));
    result_noinf = cell(size(JUMP,1),size(CATCH,1));

    #! Setup things to save
    S_H_inf 	= zeros(length(JUMP),length(CATCH));
    S_TS_inf	= zeros(length(JUMP),length(CATCH));
    S_H_noinf 	= zeros(length(JUMP),length(CATCH));
    S_TS_noinf	= zeros(length(JUMP),length(CATCH));
    S_PC_rp = zeros(length(JUMP),length(CATCH));
    S_tauh 	= zeros(length(JUMP),length(CATCH));
    S_taul 	= zeros(length(JUMP),length(CATCH));
    S_taus1 = zeros(length(JUMP),length(CATCH));

    for i = 1:length(JUMP)
        for j = 1:length(CATCH)
            print("$i $j\n")

        	##! Set params
            PRM.PS_p = JUMP[i]; # tau_l
            PRM.PC_q = CATCH[j]; # tau_h
            PRM.PC_rp = fnc_optimal_PCrp(); # optimal turning prob
            PRM.PC_lambda = 0;

            ##! Calculate taus
			tauh = PRM.PF_n/PRM.PC_q
			taul = 1./PRM.PS_p
			taus1 = fnc_taus1()
			println("tau_h/tau_s1: $(tauh/taus1) tau_l/tau_s1:  $(taul/taus1)")

			#! Save taus
            S_tauh[i,j] = tauh
            S_taul[i,j] = taul
            S_taus1[i,j] = taus1
            S_PC_rp[i,j] = PRM.PC_rp
            
            ##! Without information
            #println("Without information")
            PRM.PC_lambda=0.;
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            fnc_init_sn(cons,FLAGS);
            FLAGS["ensemble"] = 5
            FLAGS["whilecond"] = 1
            make_ensemble(school,fish,cons,fishtree,EVENTS,FLAGS);

            #save
            S_H_noinf[i,j]  = mean(cons.measure["Hrate"]);
            S_TS_noinf[i,j] = mean(cons.measure["Ts"])

            #result_noinf[i,j] = PRM.PS_p,PRM.PC_q,mean(cons.measure["Ts"]),
            #	mean(cons.measure["Hrate"]),PRM.PC_rp;

            ##! With information
            #println("With full information")
            PRM.PC_lambda=1;
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            fnc_init_sn(cons,FLAGS)
            FLAGS["ensemble"] = 5
            FLAGS["whilecond"] = 1
            make_ensemble(school,fish,cons,fishtree,EVENTS,FLAGS,3);
            result_inf[i,j] = PRM.PS_p,PRM.PC_q,mean(cons.measure["Ts"]),
            	mean(cons.measure["Hrate"]),PRM.PC_rp

            #save
            S_H_inf[i,j]  = mean(cons.measure["Hrate"]);
            S_TS_inf[i,j] = mean(cons.measure["Ts"])

        end
    end
    npzwrite("./Data/Data_2fisher.npz", ["PC_rp"=>S_PC_rp,"x"=>1,"tauh"=>S_tauh,
    									 "taul"=>S_taul,"taus1"=>S_taus1,
    									 "H_inf"=>S_H_inf,"TS_inf"=>S_TS_inf,
    									 "H_noinf"=>S_H_noinf,"TS_noinf"=>S_TS_noinf])
    #save_results(result_inf,R"PS_p PC_q \tau_s^R H PC_rp",
    #	"Data_2fisher_1school_inf",PRM)
    #save_results(result_noinf,R"PS_p PC_q \tau_s^R H PC_rp",
   # 	"Data_2fisher_1school_noinf",PRM)

end


#============== M FISHERS OPTIMAL CLIQUE SIZE ==============#
#! when we know TS_VOI is minimimzed
function do_nfisher_ts()
    #! Parameters
    PRM.PC_lambda = 1.
    PRM.PC_n = 30; # number of fishers
    
    #! Get Tau_h and Tau_l from 2fisher case
    #! find these values by plotting fig3
    #! where there is MAXIMAL VOI
    data = npzread("./Data/Data_2fisher.npz")
    taul = data["taul"]
    tauh = data["tauh"]
    TSinf = data["TS_inf"]
    TSnoinf = data["TS_noinf"]
    Z = TSinf./TSnoinf
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

    #! Iterate across different social networks
    for i = 1:300 # number of trials

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

        ## Running statistics (average Hrate within a clique)
        for p in 1:length(part)
            #for each element in the partition
            occur[part[p]]+=1
            x=(cons.measure["Hrate"][cliq[p]])
            Hrate[part[p],:]+=[mean(x) mean(x)^2]
            x=(cons.measure["Hdist"][cliq[p]])
            Hdist[part[p],:]+=[mean(x) mean(x)^2]
            x=(cons.measure["Ts"][cliq[p]])
            TS[part[p],:]+=[mean(x) mean(x)^2]
        end
    end
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
	npzwrite("./Data/Data_nfisher_ts.npz", ["Hrate_mu"=>Hrate[:,1],"x"=>1,
			"Hdist_mu"=>Hdist[:,1],"TS_mu"=>TS[:,1],"Hrate_var"=>Hrate[:,2],
            "Hdist_var"=>Hdist[:,2],"TS_var"=>TS[:,2],"occur"=>occur])

    #! Arrange results
    #coll=hcat(Hrate[:,1],Hdist[:,1],TS[:,1],Hrate[:,2],Hdist[:,2],TS[:,2],occur)
    #result=cell(size(coll,1))
    #for i =1:size(coll,1)
    #    result[i]=coll[i,:]
    #end

	##! Save
    #save_results(result,R"Hrate Hdist TS Hrate_var Hdist_var TS_var occur",
    #	"Data_rndcliq_max",PRM)
end


#============== M FISHERS OPTIMAL CLIQUE SIZE ==============#
#! when we know VOI is maximal
function do_nfisher_max()
    println("Optimal clique size vs all other sizes for 1 tauh and 1 taul")

    #! Parameters
    PRM.PC_lambda = 1.
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

    #! Iterate across different social networks
    for i = 1:300 # number of trials

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

        ## Running statistics (average Hrate within a clique)
        for p in 1:length(part)
            #for each element in the partition
            occur[part[p]]+=1
            x=(cons.measure["Hrate"][cliq[p]])
            Hrate[part[p],:]+=[mean(x) mean(x)^2]
            x=(cons.measure["Hdist"][cliq[p]])
            Hdist[part[p],:]+=[mean(x) mean(x)^2]
            x=(cons.measure["Ts"][cliq[p]])
            TS[part[p],:]+=[mean(x) mean(x)^2]
        end
    end
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
	npzwrite("./Data/Data_nfisher_max.npz", ["Hrate_mu"=>Hrate[:,1],"x"=>1,
			"Hdist_mu"=>Hdist[:,1],"TS_mu"=>TS[:,1],"Hrate_var"=>Hrate[:,2],
            "Hdist_var"=>Hdist[:,2],"TS_var"=>TS[:,2],"occur"=>occur])

    #! Arrange results
    #coll=hcat(Hrate[:,1],Hdist[:,1],TS[:,1],Hrate[:,2],Hdist[:,2],TS[:,2],occur)
    #result=cell(size(coll,1))
    #for i =1:size(coll,1)
    #    result[i]=coll[i,:]
    #end

	##! Save
    #save_results(result,R"Hrate Hdist TS Hrate_var Hdist_var TS_var occur",
    #	"Data_rndcliq_max",PRM)
end


#============== M FISHERS OPTIMAL CLIQUE SIZE ==============#
#! when we know VOI is middling
function do_nfisher_quantile(Q)
	#! Q is the quantile you want to investigate
	#! data file is saved with this number
    #! println("Optimal clique size vs all other sizes for 1 tauh and 1 taul")

    #! Parameters
    PRM.PC_lambda = 1.
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

	#! find medianish VOI
	q = quantile(Z[:],Q); #<-- Q is used here
	qmin = Z.-q;
	id = find(abs(qmin) .== minimum(abs(qmin)))[1]
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

    #! Iterate across different social networks
    for i = 1:300 # number of trials

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

        ## Running statistics (average Hrate within a clique)
        for p in 1:length(part)
            #for each element in the partition
            occur[part[p]]+=1
            x=(cons.measure["Hrate"][cliq[p]])
            Hrate[part[p],:]+=[mean(x) mean(x)^2]
            x=(cons.measure["Hdist"][cliq[p]])
            Hdist[part[p],:]+=[mean(x) mean(x)^2]
            x=(cons.measure["Ts"][cliq[p]])
            TS[part[p],:]+=[mean(x) mean(x)^2]
        end
    end
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
	ti = string(int(Q*100)+1000)
	npzwrite("./Data/Data_nfisher_"ti[2:end]".npz", ["Hrate_mu"=>Hrate[:,1],"x"=>1,
			"Hdist_mu"=>Hdist[:,1],"TS_mu"=>TS[:,1],"Hrate_var"=>Hrate[:,2],
            "Hdist_var"=>Hdist[:,2],"TS_var"=>TS[:,2],"occur"=>occur])

    #! Arrange results
    #coll=hcat(Hrate[:,1],Hdist[:,1],TS[:,1],Hrate[:,2],Hdist[:,2],TS[:,2],occur)
    #result=cell(size(coll,1))
    #for i =1:size(coll,1)
    #    result[i]=coll[i,:]
    #end

	##! Save
    #save_results(result,R"Hrate Hdist TS Hrate_var Hdist_var TS_var occur",
    #	"Data_rndcliq_max",PRM)
end



#============== M FISHERS OPTIMAL CLIQUE SIZE ==============#
#! when we know VOI is minimal
function do_nfisher_min()
    println("Optimal clique size vs all other sizes for 1 tauh and 1 taul")

    #! Parameters
    PRM.PC_lambda = 1.
    PRM.PC_n = 30; # number of fishers
    
    #! Get Tau_h and Tau_l from 2fisher case
    #! find these values by plotting fig3
    #! where there is MAXIMAL VOI
    data = npzread("./Data/Data_2fisher.npz")
    taul = data["taul"]
    tauh = data["tauh"]
    Hinf = data["H_inf"] 
    Hnoinf = data["H_noinf"]
    Z = Hinf ./ Hnoinf
    S_p = 1 ./ taul
    C_q = PRM.PF_n ./ tauh

	id = find(Z.==minimum(Z))[1]
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

    #! Iterate across different social networks
    for i = 1:300 # number of trials

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

        ## Running statistics (average Hrate within a clique)
        for p in 1:length(part)
            #for each element in the partition
            occur[part[p]]+=1
            x=(cons.measure["Hrate"][cliq[p]])
            Hrate[part[p],:]+=[mean(x) mean(x)^2]
            x=(cons.measure["Hdist"][cliq[p]])
            Hdist[part[p],:]+=[mean(x) mean(x)^2]
            x=(cons.measure["Ts"][cliq[p]])
            TS[part[p],:]+=[mean(x) mean(x)^2]
        end
    end
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
	npzwrite("./Data/Data_nfisher_min.npz", ["Hrate_mu"=>Hrate[:,1],"x"=>1,
			"Hdist_mu"=>Hdist[:,1],"TS_mu"=>TS[:,1],"Hrate_var"=>Hrate[:,2],
            "Hdist_var"=>Hdist[:,2],"TS_var"=>TS[:,2],"occur"=>occur])

    #! Arrange results
    #coll=hcat(Hrate[:,1],Hdist[:,1],TS[:,1],Hrate[:,2],Hdist[:,2],TS[:,2],occur)
    #result=cell(size(coll,1))
    #for i =1:size(coll,1)
    #    result[i]=coll[i,:]
    #end

	##! Save
    #save_results(result,R"Hrate Hdist TS Hrate_var Hdist_var TS_var occur",
    #	"Data_rndcliq_max",PRM)
end
















#============= OPTIMAL CLIQUE SIZE ================#
#! Over combinations of tau_h and tau_l
function do_Mfishers_Nschools_cliq_b()
    println("Optimal clique size vs all other sizes for every tauh taul")
    PRM.PC_lambda = 1.
    PRM.PC_n = 30;
    PRM.PS_p = 0.001;
    PRM.PC_q = 1.;
    PRM.PF_n = 1000;
    PRM.PC_f=PRM.GRD_mx*0.03
    PRM.PF_sig=PRM.GRD_mx*0.04
    PRM.PS_n = int(round(.05 / (pi *PRM.PF_sig^2 /PRM.GRD_nx^2 ) )) ;
    println("Number of schools: $(PRM.PS_n)")
    
    JUMP = logspace(log10(0.001),log10(.1),20);
    CATCH  = logspace(log10(0.01),log10(1),20);
    LAM  = linspace(0,1,10);
    allresult = cell(size(JUMP,1),size(CATCH,1));
    
    for i = 1:length(JUMP)
        for j = 1:length(CATCH)
            PRM.PS_p = JUMP[i];
            PRM.PC_q = CATCH[j];
            PRM.PC_rp = fnc_optimal_PCrp();
            tauh=PRM.PF_n/PRM.PC_q
            taul=1./PRM.PS_p
            taus=fnc_taus1()
    
            print("$i $j \n")
            println("tau_h/tau_s: $(tauh/taus) tau_l/tau_s:  $(taul/taus)")
            
            occur=zeros(PRM.PC_n) #occurrence of a clique size
            H=zeros(PRM.PC_n,2) #average catch rate per clique size
            Hdist=zeros(PRM.PC_n,2) #average CPUE per clique size
            TS=zeros(PRM.PC_n,2) #average search time per clique size
            for i = 1:300
                #random partition of PC_n
                part=pycall(pypartition["random_partition"],PyAny,PRM.PC_n) 
                cliq=cell(length(part)) #composition of the cliques
                idx=1
                for p in 1:length(part)
                    cliq[p]=idx:idx+part[p]-1
                    idx+=part[p]
                end
                #Run
                println(i)
                school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
                fnc_init_sn(cons,FLAGS,cliq)
                make_ensemble(school,fish,cons,fishtree,EVENTS,FLAGS,1);
                for p in 1:length(part)
                    #for each element in the partition
                    occur[part[p]]+=1
                    x=(cons.measure["Hrate"][cliq[p]])
                    H[part[p],:]+=[mean(x) mean(x)^2]
                    x=(cons.measure["Hdist"][cliq[p]])
                    Hdist[part[p],:]+=[mean(x) mean(x)^2]
                    x=(cons.measure["Ts"][cliq[p]])
                    TS[part[p],:]+=[mean(x) mean(x)^2]
                end
            end
            for p=1:PRM.PC_n
                if occur[p]>0
                    H[p,:]/=occur[p]
                    Hdist[p,:]/=occur[p]
                    TS[p,:]/=occur[p]
                    #Variances
                    H[p,2]-=H[p,1]^2
                    Hdist[p,2]-=Hdist[p,1]^2
                    TS[p,2]-=TS[p,1]^2
                end
            end
            coll=hcat(H[:,1],Hdist[:,1],TS[:,1],H[:,2],Hdist[:,2],TS[:,2],occur)
            result=cell(size(coll,1))
            for i =1:size(coll,1)
                result[i]=coll[i,:]
            end
        end
    end    
    
    save_results(result,
    	R"PS_p PC_q PC_rp Hrate Hdist TS Hrate_var Hdist_var TS_var occur",
    	 "Data_rndcliq_b",PRM)
end






function do_Mfishers_Nschools_cliq_b()
    println("Optimal clique size vs all other sizes for every tauh taul")
    PRM.PC_lambda = 1.
    PRM.PC_n = 30;
    PRM.PS_p = 0.001;
    PRM.PC_q = 1.;
    PRM.PF_n = 1000;
    PRM.PC_f=PRM.GRD_mx*0.03
    PRM.PF_sig=PRM.GRD_mx*0.04
    PRM.PS_n = int(round(.05 / (pi *PRM.PF_sig^2 /PRM.GRD_nx^2 ) )) ;
    println("Number of schools: $(PRM.PS_n)")
    
    JUMP = logspace(log10(0.001),log10(.1),20);
    CATCH  = logspace(log10(0.01),log10(1),20);
    LAM  = linspace(0,1,10);
    allresult = cell(size(JUMP,1),size(CATCH,1));
    
    for i = 1:length(JUMP)
        for j = 1:length(CATCH)
            PRM.PS_p = JUMP[i];
            PRM.PC_q = CATCH[j];
            PRM.PC_rp = fnc_optimal_PCrp();
            tauh=PRM.PF_n/PRM.PC_q
            taul=1./PRM.PS_p
            taus=fnc_taus1()
    
            print("$i $j \n")
            println("tau_h/tau_s: $(tauh/taus) tau_l/tau_s:  $(taul/taus)")
            
            occur=zeros(PRM.PC_n) #occurrence of a clique size
            H=zeros(PRM.PC_n,2) #average catch rate per clique size
            Hdist=zeros(PRM.PC_n,2) #average CPUE per clique size
            TS=zeros(PRM.PC_n,2) #average search time per clique size
            for i = 1:300
                #random partition of PC_n
                part=pycall(pypartition["random_partition"],PyAny,PRM.PC_n) 
                cliq=cell(length(part)) #composition of the cliques
                idx=1
                for p in 1:length(part)
                    cliq[p]=idx:idx+part[p]-1
                    idx+=part[p]
                end
                #Run
                println(i)
                school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
                fnc_init_sn(cons,FLAGS,cliq)
                make_ensemble(school,fish,cons,fishtree,EVENTS,FLAGS,1);
                for p in 1:length(part)
                    #for each element in the partition
                    occur[part[p]]+=1
                    x=(cons.measure["Hrate"][cliq[p]])
                    H[part[p],:]+=[mean(x) mean(x)^2]
                    x=(cons.measure["Hdist"][cliq[p]])
                    Hdist[part[p],:]+=[mean(x) mean(x)^2]
                    x=(cons.measure["Ts"][cliq[p]])
                    TS[part[p],:]+=[mean(x) mean(x)^2]
                end
            end
            for p=1:PRM.PC_n
                if occur[p]>0
                    H[p,:]/=occur[p]
                    Hdist[p,:]/=occur[p]
                    TS[p,:]/=occur[p]
                    #Variances
                    H[p,2]-=H[p,1]^2
                    Hdist[p,2]-=Hdist[p,1]^2
                    TS[p,2]-=TS[p,1]^2
                end
            end
            coll=hcat(H[:,1],Hdist[:,1],TS[:,1],H[:,2],Hdist[:,2],TS[:,2],occur)
            result=cell(size(coll,1))
            for i =1:size(coll,1)
                result[i]=coll[i,:]
            end
        end
    end    
    
    save_results(result,
    	R"PS_p PC_q PC_rp Hrate Hdist TS Hrate_var Hdist_var TS_var occur",
    	 "Data_rndcliq_b",PRM)
end



