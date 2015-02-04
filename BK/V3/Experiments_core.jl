
##########################!!!!! SIMULATION EXPERIMENTS !!!!!!!!################


##################### BASIC TESTS
#! Timing test: ... this does what
function do_timingtest(fsave=true)
    println("Timing test")

    PRM.PC_lambda = 1.
    PRM.PC_n = 10;
    PRM.PC_f=PRM.GRD_mx*0.03
    PRM.PF_sig=PRM.GRD_mx*0.04
    PRM.PS_p = 0.001;
    PRM.PC_q = 1.;
    PRM.PS_n = 10;
    PRM.PF_n = 100;
    
    PRM.PC_rp = fnc_optimal_PCrp();
    tauh=PRM.PF_n/PRM.PC_q
    taul=1./PRM.PS_p
    taus=fnc_taus1()
    println("tau_h/tau_s: $(tauh/taus) tau_l/tau_s:  $(taul/taus)")
    part=pycall(pypartition["random_partition"],PyAny,PRM.PC_n)
    cliq=cell(length(part)) #composition of the cliques
    idx=1
    for p in 1:length(part)
        cliq[p]=idx:idx+part[p]-1
        idx+=part[p]
    end

    school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
    fnc_init_sn(cons,FLAGS,cliq)
    FLAGS["save"]=fsave
    @time  make_season(school,fish,cons,fishtree,EVENTS,FLAGS,1,OUT);
end


##################### 1 SCHOOl 1 FISHER
function do_1fisher_1school_a()
    ######  Expected search time (T_s^*)
    ##! designed to link up to Benichou analytics
    println("Expected search time")
    RP = linspace(0.2,.9,20); #! probability of turning per unit time
    NN = 10; #! number of replicates
    PRM.PC_n=1 #! number of fishers
    result = cell(length(RP),NN); #! array for results collection

    for i = 1:length(RP) # Iterate over different turning probs
    	for j = 1:NN # run sim a number of times 
			PRM.PC_rp = RP[i];
			school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
			FLAGS["measure_frac"]=true
			fnc_init_sn(cons,FLAGS)
			make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
			result[i,j] = RP[i],mean(cons.measure["Ts"]),
						mean(cons.measure["Tv"]),mean(cons.measure["f1"]);
			print(i/length(RP)," | ", j/NN,"\n")
		end
    end
    save_results(result,R"PC_rp \tau_s^R std f1","Data_1fisher_1school_a",PRM)
end

function do_1fisher_1school_b()
    ###### T_s^* as a function of C_f and F_sig
    SIG = linspace(2,GRD_mx/10,32);
    FF  = linspace(2,GRD_mx/10,32);
    result = cell(size(SIG,1),size(FF,1));
    PRM.PC_n=1
    for i = 1:length(SIG)
        for j = 1:length(FF)
            PRM.PF_sig = SIG[i];
            PRM.PC_f   = FF[j];
            PRM.PC_rp = fnc_optimal_PCrp();
            print(i," ",j,"\n")
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            fnc_init_sn(cons,FLAGS)
            FLAGS["measure_frac"]=true
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            result[i,j] =SIG[i],FF[j], mean(cons.measure["Ts"]), 
            	mean(cons.measure["Tv"]),mean(cons.measure["f1"]),PRM.PC_rp;
        end
    end
    save_results(result,R"PF_sig PC_f \tau_s^R std f1 PC_rp",
    	"Data_1fisher_1school_b",PRM)
end


##################### 1 SCHOOL 2 FISHERS
##! Test effect of friendship, on Tau_s and CPUE, 
##! as Tau_l and Tau_h are varied
function do_2fisher_1school()
    ###### Test performance as a function of tau_l and tau_h
    ## basic values of C_f and F_sig
    PRM.PC_f=PRM.GRD_mx*0.03
    PRM.PF_sig=PRM.GRD_mx*0.04
    PRM.PC_n = 2; #2 fishers!
    PRM.PS_n = int(round(.05 / (pi *PRM.PF_sig^2 /PRM.GRD_nx^2 ) )) ;
    println("Performance as a function of tau_l and tau_h")
    println("Number of schools: $(PRM.PS_n)")
    
    #! Choose resolution of Tau_l and Tau_h
    JUMP = logspace(log10(0.001),log10(.1),20); ## Tau_l
    CATCH  = logspace(log10(0.01),log10(1),20); ## Tau_h
    result_inf = cell(size(JUMP,1),size(CATCH,1));
    result_noinf = cell(size(JUMP,1),size(CATCH,1));

    for i = 1:length(JUMP)
        for j = 1:length(CATCH)
            PRM.PS_p = JUMP[i]; # tau_l
            PRM.PC_q = CATCH[j]; # tau_h
            PRM.PC_rp = fnc_optimal_PCrp(); # optimal turning prob
            print("$i $j\n")
            
            ##! Without information
            println("Without information")
            PRM.PC_lambda=0;
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            fnc_init_sn(cons,FLAGS)
            FLAGS["measure_frac"]=true
            FLAGS["measure_fij"]=true
            
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            result_noinf[i,j] = PRM.PS_p,PRM.PC_q,mean(cons.measure["Ts"]),
            	mean(cons.measure["f1"]),mean(cons.measure["f2"]),
            	mean(cons.measure["fij"]),mean(cons.measure["bound"]),
            	mean(cons.measure["Hrate"]),PRM.PC_rp;

            ##! With information
            println("With full information")
            PRM.PC_lambda=1;
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            fnc_init_sn(cons,FLAGS)
            FLAGS["measure_frac"]=true
            FLAGS["measure_fij"]=true
            
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            result_inf[i,j] = PRM.PS_p,PRM.PC_q,mean(cons.measure["Ts"]),
            	mean(cons.measure["f1"]),mean(cons.measure["f2"]),
            	mean(cons.measure["fij"]),mean(cons.measure["bound"]),
            	mean(cons.measure["Hrate"]),PRM.PC_rp

        end
    end
    save_results(result_inf,R"PS_p PC_q \tau_s^R f1 f2 fij bound H PC_rp",
    	"Data_2fisher_1school_inf",PRM)
    save_results(result_noinf,R"PS_p PC_q \tau_s^R f1 f2 fij bound H PC_rp",
    	"Data_2fisher_1school_noinf",PRM)

end


##################### M FISHERS N SCHOOL 
##! Choose M Fishers such that different clique sizes are possible
##! Choose N Schools such that the area of the domain covered by fish
##! at any one moment is ~ 5%
function do_mfisher_nschool_cliq_a()
    println("Optimal clique size vs all other sizes for 1 tauh and 1 taul")
    PRM.PC_lambda = 1.
    PRM.PC_n = 30;
    PRM.PC_f=PRM.GRD_mx*0.04
    PRM.PF_sig=PRM.GRD_mx*0.04
    PRM.PS_p = 0.001;
    PRM.PC_q = 1.;
    PRM.PS_n = 10;
    PRM.PF_n = 1000;
    
    PRM.PC_rp = fnc_optimal_PCrp();
    tauh=PRM.PF_n/PRM.PC_q
    taul=1./PRM.PS_p
    taus=fnc_taus1()
    println("tau_h/tau_s: $(tauh/taus) tau_l/tau_s:  $(taul/taus)")
    
    occur=zeros(PRM.PC_n) #occurrence of a clique size
    H=zeros(PRM.PC_n,2) #average catch rate per clique size
    Hdist=zeros(PRM.PC_n,2) #average CPUE per clique size
    TS=zeros(PRM.PC_n,2) #average search time per clique size
    for i = 1:300 # number of trials
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
        make_season(school,fish,cons,fishtree,EVENTS,FLAGS,1);
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
    save_results(result,R"Hrate Hdist TS Hrate_var Hdist_var TS_var occur",
    	"Data_rndcliq",PRM)
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
                make_season(school,fish,cons,fishtree,EVENTS,FLAGS,1);
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
    	 "Data_rndcliq_explor",PRM)
end



