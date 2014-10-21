
##########################!!!!! SIMULATION EXPERIMENTS !!!!!!!!################


##################### BASIC TESTS
    
function do_timingtest()
    println("Timing test")
    PRM.PC_rp = 0.72; # choose random change in walk
    school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
    FLAGS["save"]=true
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
    for i = 1:length(RP)
        PRM.PC_rp = RP[i];
        time=[]
        for iter in 1:300
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
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
    RP = 1:2:60 #[1, 2,4,9,16,25,36,49,64];
    Ts = cell(size(RP));
    PRM.PS_p=0.0
    for i = 1:length(RP)
        PRM.PS_n = RP[i];
        time=[]
        dist=[]
        for iter in 1:300
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
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
    result = cell(size(RP));

    for i = 1:length(RP)
        PRM.PC_rp = RP[i];
        school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
        FLAGS["measure_frac"]=true
        make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
        result[i] = RP[i],mean(cons.measure["Ts"]), mean(cons.measure["Tv"]),mean(cons.measure["f1"]);
        print(i,"\n")
    end
    save_results(result,R"PC_rp \tau_s^R std f1","Data_Fig2a",PRM)
end

function do_fig2b()
    ###### Test performance as a function of C_f and F_sig
    SIG = linspace(2,GRD_mx/10,32);
    FF  = linspace(2,GRD_mx/10,32);
    result = cell(size(SIG,1),size(FF,1));
    for i = 1:length(SIG)
        for j = 1:length(FF)
            PRM.PF_sig = SIG[i];
            PRM.PC_f   = FF[j];
            PRM.PC_rp = fnc_optimal_PCrp();
            print(i," ",j,"\n")
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            FLAGS["measure_frac"]=true
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            result[i,j] =SIG[i],FF[j], mean(cons.measure["Ts"]), mean(cons.measure["Tv"]),mean(cons.measure["f1"]),PRM.PC_rp;
        end
    end
    save_results(result,R"PF_sig PC_f \tau_s^R std f1 PC_rp","Data_Fig2b",PRM)
end


##################### 1 SCHOOL 2 FISHERS
######## Test effect of friendship, on Tau_s and CPUE, as Tau_l is varied

function do_fig3()
    ###### Test performance as a function of tau_l and tau_h
    ## basic values of C_f and F_sig
    PRM.PC_f=PRM.GRD_nx*0.03
    PRM.PF_sig=PRM.GRD_nx*0.04
    PRM.PC_n = 2; #2 fishers!
    PRM.PS_n = int(round(.10 / (pi *PRM.PF_sig^2 /PRM.GRD_nx^2 ) )) ;
    println("Performance as a function of tau_l and tau_h")
    println("Number of schools: $(PRM.PS_n)")
    
    JUMP = logspace(log10(0.001),log10(.1),32);
    CATCH  = logspace(log10(0.01),log10(1),32);
    result_inf = cell(size(JUMP,1),size(CATCH,1));
    result_noinf = cell(size(JUMP,1),size(CATCH,1));

    for i = 1:length(JUMP)
        for j = 1:length(CATCH)
            PRM.PS_p = JUMP[i];
            PRM.PC_q = CATCH[j];
            PRM.PC_rp = fnc_optimal_PCrp();
            print("$i $j \n")
            
            println("Without information")
            #Without information
            PRM.PC_lambda=0;
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            FLAGS["measure_frac"]=true
            
#            println(cons.SN)
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            result_noinf[i,j] =PRM.PS_p,PRM.PC_q,mean(cons.measure["Ts"]),mean(cons.measure["f1"]),mean(cons.measure["fij"]),mean(cons.measure["Hrate"]);

            println("With full information")
            #With information
            PRM.PC_lambda=1;
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            FLAGS["measure_frac"]=true
            
#            println(cons.SN)
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);

            result_inf[i,j] =PRM.PS_p,PRM.PC_q,mean(cons.measure["Ts"]),mean(cons.measure["f1"]),mean(cons.measure["fij"]),mean(cons.measure["Hrate"]),PRM.PC_rp

        end
    end
    save_results(result_inf,R"PS_p PC_q \tau_s^R f1 fij H PC_rp","Data_Fig3_inf",PRM)
    save_results(result_noinf,R"PS_p PC_q \tau_s^R f1 fij H PC_rp","Data_Fig3_noinf",PRM)

end


function do_fig4opt()
    println("Optimal lambda")
    RP = linspace(0.,1.,5);
    NS = [1,2,3,5,10,20]; #Number of schools
    result = cell(size(RP,1),size(NS,1));

    PRM.PC_n = 10;
    PRM.PC_ncliq = 1;
    PRM.PC_f=PRM.GRD_nx*0.05
    PRM.PF_sig=PRM.GRD_nx*0.05
    PRM.PS_p = 0.001;
    PRM.PC_q = 0.01;
    PRM.PC_rp = fnc_optimal_PCrp();
    for i = 1:length(RP)
        for j=1:length(NS)
            PRM.PC_lambda = RP[i];
            PRM.PS_n = NS[j];
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            FLAGS["measure_frac"]=true
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            result[i,j] = RP[i], NS[j],mean(cons.measure["Ts"]), mean(cons.measure["Hrate"]),mean( cons.measure["Hdist"]),mean( cons.measure["f2"]);
            println("$i $j")
        end
    end
    save_results(result,R"PC_lambda PS_n \tau_s^R Hrate Hdist f2","Data_Fig4opt",PRM)

end

function do_fig4opt_cliq()
    println("Optimal number of cliques")
    RP = [1,2,3,5,10];
    NS = [1,2,3,5,10,20]; #Number of schools
    result = cell(size(RP,1),size(NS,1));

    PRM.PC_n = 10;
    PRM.PC_f=PRM.GRD_nx*0.05
    PRM.PF_sig=PRM.GRD_nx*0.05
    PRM.PS_p = 0.001;
    PRM.PC_q = 0.01;
    PRM.PF_n = 50;
    PRM.PC_lambda = 1.;
    PRM.PC_rp = fnc_optimal_PCrp();
    for i = 1:length(RP)
        for j=1:length(NS)
            PRM.PC_ncliq = RP[i];
            PRM.PS_n = NS[j];
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            FLAGS["measure_frac"]=true
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            result[i,j] = RP[i], NS[j],mean(cons.measure["Ts"]),mean(cons.measure["Hrate"]),mean( cons.measure["Hdist"]),mean( cons.measure["f2"]);
            println("$i $j")
        end
    end
    save_results(result,R"PC_ncliq PS_n \tau_s^R Hrate Hdist f2","Data_Fig4opt_cliq",PRM)
end


function do_rndcliq()
    println("Optimal clique size vs all other clique sizes")
    PRM.PC_n=30
    PRM.PC_f=PRM.GRD_nx*0.04
    PRM.PF_sig=PRM.GRD_nx*0.04
    PRM.PS_n=20
    PRM.PC_lambda = 1.

    occur=zeros(PRM.PC_n) #occurrence of a clique size
    H=zeros(PRM.PC_n) #average catch rate per clique size
    TS=zeros(PRM.PC_n) #average search time per clique size
    for i = 1:1000
        part=pycall(pypartition["random_partition"],PRM.PC_n) #random partition of PC_n
        cliq=cell(length(part)) #composition of the cliques
        idx=1
        for p in part
            cliq[p]=idx:idx+p-1
            idx+=p
        end
        #Run
        school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium(cliques=cliq);
        FLAGS["measure_frac"]=true
        make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
        for p in part
            #for each element in the partition
            occur[p]+=1
            H[p]+=mean(cons.measure["Hrate"][cliq[p]])
            TS[p]+=mean(cons.measure["TS"][cliq[p]])
        end
    end
    H/=occur
    TS/=occur
    result=hcat(1:PC_n,H,TS,occur)'
    save_results(result,R"size H TS occur","Data_Fig4opt_cliq",PRM)
end


######## Test effect of friendship, on Tau_s and CPUE, as F_n is varied

##################### 1 FISH N FISHERS (do 5, 10, 20 fishers)
######## Stupid Optimal Social Network test (just keep adding links)

######## Evolving social network - find the ESS 








######## END  #######


