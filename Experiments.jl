
##########################!!!!! SIMULATION EXPERIMENTS !!!!!!!!################


##################### BASIC TESTS
    
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
    part=pycall(pypartition["random_partition"],PyAny,PRM.PC_n) #random partition of PC_n
    cliq=cell(length(part)) #composition of the cliques
    idx=1
    for p in 1:length(part)
        cliq[p]=idx:idx+part[p]-1
        idx+=part[p]
    end

    school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
    init_network(cons,FLAGS,cliq)
    FLAGS["save"]=fsave
    FLAGS["measure_frac"]=true
    @time  make_season(school,fish,cons,fishtree,EVENTS,FLAGS,1,OUT);
    #npzwrite("./Data/Data_fish.npy", OUT.fish_xy)
    #npzwrite("./Data/Data_fishers.npy", OUT.cons_xy)
    #npzwrite("./Data/Data_clusters.npy", OUT.schl_xy)
    #npzwrite("./Data/Data_harvest.npy", OUT.cons_H)
end


##################### 1 SCHOOl 1 FISHER
    
function do_first_passage()
    ###### FIRST PASSAGE time (stop as soon as a school is found)
    println("First passage time")
    RP = linspace(0.1,.95,20);
    Ts = cell(size(RP));
    PRM.PS_p=0.0
    PRM.PC_n=1
    for i = 1:length(RP)
        PRM.PC_rp = RP[i];
        time=[]
        for iter in 1:300
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            init_network(cons,FLAGS)
            time=vcat(time, make_season(school,fish,cons,fishtree,EVENTS,FLAGS,4) ) ;
        end
        Ts[i] = RP[i],mean(time),std(time);#cons.measure["Ts"], cons.measure["Tv"];
        print(i,"\n")
    end
    save_results(Ts,R"PC_rp \tau_s^R std","Data_firstpass",PRM)
end

function do_fstpass_nschool()
    ###### FIRST PASSAGE time (stop as soon as a school is found)
    ### As a function of number of schools
    println("First passage time")
    PRM.PC_n=1
    RP = 1:2:60 #[1, 2,4,9,16,25,36,49,64];
    Ts = cell(size(RP));
    PRM.PS_p=0.0
    for i = 1:length(RP)
        PRM.PS_n = RP[i];
        time=[]
        dist=[]
        for iter in 1:300
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            init_network(cons,FLAGS)
            #pos=cons.x[1,:]
            time=vcat(time, make_season(school,fish,cons,fishtree,EVENTS,FLAGS,4) ) ;
            dist=vcat(dist, cons.measure["distance"][1])#fnc_dist(cons.x[1,:],pos))
        end
        Ts[i] = RP[i], mean(time),std(time),mean(dist);#cons.measure["Ts"], cons.measure["Tv"];
        print("$i $(RP[i]) schools ($(3.14159*PRM.PF_sig^2 * RP[i]/GRD_mx^2 *100)\% covering) \n")
    end
    save_results(Ts,R"PS_n \tau_s^R std, dist","Data_firstpass_ns",PRM)
end

function do_fig2a()
    ######  TYPICAL SEARCH time over many searches
    println("Typical search time")
    RP = linspace(0.2,.9,20);
    PRM.PC_n=1
    result = cell(size(RP));

    for i = 1:length(RP)
        PRM.PC_rp = RP[i];
        school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
        FLAGS["measure_frac"]=true
        init_network(cons,FLAGS)
        make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
        result[i] = RP[i],mean(cons.measure["Ts"]), mean(cons.measure["Tv"]),mean(cons.measure["f1"]);
        print(i,"\n")
    end
    save_results(result,R"PC_rp \tau_s^R std f1","Data_Fig2a",PRM)
end


function do_spying()
    ######  TYPICAL SEARCH Time as a function of spying radius
    println("Everybody is spying")
    SPY = linspace(PRM.PC_f,PRM.GRD_mx/2,20);
    PRM.PC_n=10
    PRM.PC_lambda=1.
    result = cell(size(SPY));
    
    println("Movie")
    school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
    PRM.PC_spy=10
    FLAGS["spying"]=true
    FLAGS["save"]=true
    init_network(cons,FLAGS)
    make_season(school,fish,cons,fishtree,EVENTS,FLAGS,1,OUT);

    
    println("Main run")
    for i = 1:length(SPY)
        PRM.PC_spy = SPY[i];
        school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
        FLAGS["spying"]=true
        init_network(cons,FLAGS)
        make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
        result[i] = SPY[i],mean(cons.measure["Ts"]), mean(cons.measure["Hrate"][1:PRM.PC_n/2]), mean(cons.measure["Hrate"][PRM.PC_n/2:end]);
        print(i,"\n")
    end
    save_results(result,R"PC_spy tau_s H1 H2","Data_spy",PRM)
end


function do_fig2b()
    ###### Test performance as a function of C_f and F_sig
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
            init_network(cons,FLAGS)
            FLAGS["measure_frac"]=true
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            result[i,j] =SIG[i],FF[j], mean(cons.measure["Ts"]), mean(cons.measure["Tv"]),mean(cons.measure["f1"]),PRM.PC_rp;
        end
    end
    save_results(result,R"PF_sig PC_f \tau_s^R std f1 PC_rp","Data_Fig2b",PRM)
end

function do_fig2c()
    ###### Test performance as a function of jump and catch proba
    PRM.PC_f=PRM.GRD_mx*0.03
    PRM.PF_sig=PRM.GRD_mx*0.04
    PRM.PS_n = 10 ;
    PRM.PF_n=10
    PRM.PC_n=1
    JUMP = logspace(log10(0.001),log10(.1),32);
    CATCH  = logspace(log10(0.01),log10(1),32);
    result = cell(size(JUMP,1),size(CATCH,1));
    
    for i = 1:length(JUMP)
        for j = 1:length(CATCH)
            PRM.PS_p = JUMP[i];
            PRM.PC_q = CATCH[j];
            PRM.PC_rp = fnc_optimal_PCrp();
            print("$i $j\n")
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            init_network(cons,FLAGS)
            FLAGS["measure_frac"]=true
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            result[i,j] =JUMP[i],CATCH[j], mean(cons.measure["Ts"]), mean(cons.measure["Tv"]),mean(cons.measure["f1"]),PRM.PC_rp;
        end
    end
    save_results(result,R"PS_p PC_q \tau_s^R std f1 PC_rp","Data_Fig2c",PRM)
end


##################### 1 SCHOOL 2 FISHERS
######## Test effect of friendship, on Tau_s and CPUE, as Tau_l is varied

function do_fig3()
    ###### Test performance as a function of tau_l and tau_h
    ## basic values of C_f and F_sig
    PRM.PC_f=PRM.GRD_mx*0.03
    PRM.PF_sig=PRM.GRD_mx*0.04
    PRM.PC_n = 2; #2 fishers!
    PRM.PS_n = int(round(.05 / (pi *PRM.PF_sig^2 /PRM.GRD_nx^2 ) )) ;
    println("Performance as a function of tau_l and tau_h")
    println("Number of schools: $(PRM.PS_n)")
    
    JUMP = logspace(log10(0.001),log10(.1),8);
    CATCH  = logspace(log10(0.01),log10(1),8);
    result_inf = cell(size(JUMP,1),size(CATCH,1));
    result_noinf = cell(size(JUMP,1),size(CATCH,1));

    for i = 1:length(JUMP)
        for j = 1:length(CATCH)
            PRM.PS_p = JUMP[i];
            PRM.PC_q = CATCH[j];
            PRM.PC_rp = fnc_optimal_PCrp();
            print("$i $j\n")
            
            println("Without information")
            #Without information
            PRM.PC_lambda=0;
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            init_network(cons,FLAGS)
            FLAGS["measure_frac"]=true
            FLAGS["measure_fij"]=true
            
#            println(cons.SN)
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            result_noinf[i,j] =PRM.PS_p,PRM.PC_q,mean(cons.measure["Ts"]),mean(cons.measure["f1"]),mean(cons.measure["f2"]),mean(cons.measure["fij"]),mean(cons.measure["bound"]),mean(cons.measure["Hrate"]),PRM.PC_rp;

            println("With full information")
            #With information
            PRM.PC_lambda=1;
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            init_network(cons,FLAGS)
            FLAGS["measure_frac"]=true
            FLAGS["measure_fij"]=true
            
#            println(cons.SN)
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);

            result_inf[i,j] =PRM.PS_p,PRM.PC_q,mean(cons.measure["Ts"]),mean(cons.measure["f1"]),mean(cons.measure["f2"]),mean(cons.measure["fij"]),mean(cons.measure["bound"]),mean(cons.measure["Hrate"]),PRM.PC_rp

        end
    end
    save_results(result_inf,R"PS_p PC_q \tau_s^R f1 f2 fij bound H PC_rp","Data_Fig3_inf",PRM)
    save_results(result_noinf,R"PS_p PC_q \tau_s^R f1 f2 fij bound H PC_rp","Data_Fig3_noinf",PRM)

end


function do_worst()
    ###### Find worst value of comm as a function of tau_l and tau_h
    PRM.PC_f=PRM.GRD_mx*0.03
    PRM.PF_sig=PRM.GRD_mx*0.04
    PRM.PC_n = 10 #10 fishers
    PRM.PS_n = int(round(.05 / (pi *PRM.PF_sig^2 /PRM.GRD_nx^2 ) )) ;
    println("Performance as a function of tau_l and tau_h")
    println("Number of schools: $(PRM.PS_n)")
    
    JUMP = logspace(log10(0.001),log10(.1),20);
    CATCH  = logspace(log10(0.01),log10(1),20);
    LAM  = linspace(0,1,10);
    result = cell(size(JUMP,1),size(CATCH,1),size(LAM,1));
    
    for i = 1:length(JUMP)
        for j = 1:length(CATCH)
            PRM.PS_p = JUMP[i];
            PRM.PC_q = CATCH[j];
            PRM.PC_rp = fnc_optimal_PCrp();
            print("$i $j \n")
            for k in 1:length(LAM)
                #Without information
                PRM.PC_lambda=LAM[k];
                school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
                init_network(cons,FLAGS)
                FLAGS["measure_frac"]=true
                
    #            println(cons.SN)
                make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
                result[i,j,k] =JUMP[i],CATCH[j],LAM[k],mean(cons.measure["Ts"]),mean(cons.measure["f1"]),mean(cons.measure["Hrate"]),mean(cons.measure["Hdist"]);
            end
        end
    end
    save_results(result,R"PS_p PC_q PC_lambda \tau_s^R f1 Hrate Hdist PC_rp","Data_worst",PRM)

end



function do_fig4opt()
    println("Optimal lambda")
    RP = linspace(0.,1.,5);
    NS = [1,2,3,5,10,20]; #Number of schools
    result = cell(size(RP,1),size(NS,1));

    PRM.PC_n = 2;
    PRM.PC_ncliq = 1;
    PRM.PC_f=PRM.GRD_mx*0.05
    PRM.PF_sig=PRM.GRD_mx*0.05
    PRM.PS_p = 0.001;
    PRM.PC_q = 0.1;
    PRM.PC_rp = fnc_optimal_PCrp();
    for i = 1:length(RP)
        for j=1:length(NS)
            PRM.PC_lambda = RP[i];
            PRM.PS_n = NS[j];
            taus1=fnc_taus1()
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            init_network(cons,FLAGS)
            FLAGS["measure_frac"]=true
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            result[i,j] = RP[i], NS[j],taus1, mean(cons.measure["Ts"]), mean(cons.measure["Hrate"]),mean(cons.measure["Hdist"]),mean(cons.measure["f1"]),mean(cons.measure["f2"]);
            println("$i $j")
        end
    end
    save_results(result,R"PC_lambda PS_n taus1 \tau_s^R Hrate Hdist f1 f2","Data_Fig4opt",PRM)

end


function do_fig4opt_cn()
    println("Explore Lambda & Cn")
    RP = linspace(0.,1.,10);
    NS = [x for x in 1:10] #[1,2,3,5,10]; #Number of fishers
    result = cell(size(RP,1),size(NS,1));

    PRM.PC_n = 2;
    PRM.PC_ncliq = 1;
    PRM.PC_f=PRM.GRD_mx*0.05
    PRM.PF_sig=PRM.GRD_mx*0.05
    PRM.PS_p = 0.001;
    PRM.PC_q = 0.1;
    PRM.PS_n=10;
    PRM.PC_rp = fnc_optimal_PCrp();
    for i = 1:length(RP)
        for j=1:length(NS)
            PRM.PC_lambda = RP[i];
            PRM.PC_n = NS[j];
            taus1=fnc_taus1()
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            init_network(cons,FLAGS)
            FLAGS["measure_frac"]=true
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            result[i,j] = RP[i], NS[j],taus1, mean(cons.measure["Ts"]), mean(cons.measure["Hrate"]),mean(cons.measure["Hdist"]),mean(cons.measure["f1"]),mean(cons.measure["f2"]);
            println("$i $j")
        end
    end
    save_results(result,R"PC_lambda PC_n taus1 \tau_s^R Hrate Hdist f1 f2","Data_Fig4opt_cn",PRM)

end

function do_fig4opt_cliq()
    println("Optimal number of cliques")
    RP = [1,2,3,5,10];
    NS = [1,2,3,5,10,20]; #Number of schools
    result = cell(size(RP,1),size(NS,1));

    PRM.PC_n = 10;
    PRM.PC_f=PRM.GRD_mx*0.05
    PRM.PF_sig=PRM.GRD_mx*0.05
    PRM.PS_p = 0.001;
    PRM.PC_q = 0.1;
    PRM.PF_n = 50;
    PRM.PC_lambda = 1.;
    PRM.PC_rp = fnc_optimal_PCrp();
    for i = 1:length(RP)
        for j=1:length(NS)
            PRM.PC_ncliq = RP[i];
            PRM.PS_n = NS[j];
            taus1=fnc_taus1()
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            init_network(cons,FLAGS)
            FLAGS["measure_frac"]=true
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            result[i,j] = RP[i], NS[j],taus1,mean(cons.measure["Ts"]),mean(cons.measure["Hrate"]),mean( cons.measure["Hdist"]),mean( cons.measure["f2"]);
            println("$i $j")
        end
    end
    save_results(result,R"PC_ncliq PS_n taus1 \tau_s^R Hrate Hdist f2","Data_Fig4opt_cliq",PRM)
end


function do_fig4opt_comp()
    println("Cliques vs Lambda")
    RP = [1,2,3,5,10];
    LAM = linspace(0.,1.,5);
    result = cell(size(RP,1),size(LAM,1));

    PRM.PC_n = 10;
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
    for i = 1:length(RP)
        for j=1:length(LAM)
            PRM.PC_ncliq = RP[i];
            PRM.PC_lambda = LAM[j];
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            init_network(cons,FLAGS)
            FLAGS["measure_frac"]=true
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            result[i,j] = RP[i], LAM[j],mean(cons.measure["Ts"]),mean(cons.measure["Hrate"]),mean( cons.measure["Hdist"]),mean( cons.measure["f2"]);
            println("$i $j")
        end
    end
    save_results(result,R"PC_ncliq PC_lambda \tau_s^R Hrate Hdist f2","Data_Fig4opt_comp",PRM)
end

function do_rndcliq()
    println("Optimal clique size vs all other sizes")
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
    Hdist=zeros(PRM.PC_n,2) #average catch/traveled distance per clique size
    TS=zeros(PRM.PC_n,2) #average search time per clique size
    for i = 1:300
        part=pycall(pypartition["random_partition"],PyAny,PRM.PC_n) #random partition of PC_n
        cliq=cell(length(part)) #composition of the cliques
        idx=1
        for p in 1:length(part)
            cliq[p]=idx:idx+part[p]-1
            idx+=part[p]
        end
        #Run
        println(i)
        school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
        init_network(cons,FLAGS,cliq)
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
    save_results(result,R"Hrate Hdist TS Hrate_var Hdist_var TS_var occur","Data_rndcliq",PRM)
end

function do_rndcliq_explor()
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
            Hdist=zeros(PRM.PC_n,2) #average catch/traveled distance per clique size
            TS=zeros(PRM.PC_n,2) #average search time per clique size
            for i = 1:300
                part=pycall(pypartition["random_partition"],PyAny,PRM.PC_n) #random partition of PC_n
                cliq=cell(length(part)) #composition of the cliques
                idx=1
                for p in 1:length(part)
                    cliq[p]=idx:idx+part[p]-1
                    idx+=part[p]
                end
                #Run
                println(i)
                school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
                init_network(cons,FLAGS,cliq)
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
    
    save_results(result,R"PS_p PC_q PC_rp Hrate Hdist TS Hrate_var Hdist_var TS_var occur","Data_rndcliq_explor",PRM)
end




######## Test effect of friendship, on Tau_s and CPUE, as F_n is varied

##################### 1 FISH N FISHERS (do 5, 10, 20 fishers)
######## Stupid Optimal Social Network test (just keep adding links)

######## Evolving social network - find the ESS 








######## END  #######


