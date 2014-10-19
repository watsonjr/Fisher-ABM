
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
    RP = [1, 2,4,9,16,25];
    Ts = cell(size(RP));
    PRM.PS_p=0.0
    for i = 1:length(RP)
        PRM.PS_n = RP[i];
        time=[]
        for iter in 1:300
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            time=vcat(time, make_season(school,fish,cons,fishtree,EVENTS,FLAGS,4) ) ;
        end
        Ts[i] = RP[i], mean(time),std(time);#cons.measure["Ts"], cons.measure["Tv"];
        print(i,"\n")
    end
    save_results(Ts,R"PS_n \tau_s^R std","Data_firstpass_ns",PRM)
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
            result[i,j] =SIG[i],FF[j], mean(cons.measure["Ts"]), mean(cons.measure["Tv"]),mean(cons.measure["f1"]);
        end
    end
    save_results(result,R"PF_sig PC_f \tau_s^R std f1","Data_Fig2b",PRM)
    #=
    TS = zeros(Float64,size(SIG,1),size(FF,1),PRM.PC_n)
    F1 = zeros(Float64,size(SIG,1),size(FF,1),PRM.PC_n)
    xs = zeros(Float64,size(SIG,1),size(FF,1),2)
    for i = 1:size(result,1)
        for j=1:size(result,2)
            TS[i,j,:] = result[i,j][2]
            xs[i,j,:]=result[i,j][1]
            F1[i,j,:] = result[i,j][3]
        end
    end
    npzwrite("./Data/Data_Fig2b.npy", TS)
    npzwrite("./Data/Data_Fig2b_f1.npy", TS)
    npzwrite("./Data/Data_Fig2b_xs.npy", xs)
    =#
end


##################### 1 SCHOOL 2 FISHERS
######## Test effect of friendship, on Tau_s and CPUE, as Tau_l is varied

function do_fig3()
    ###### Test performance as a function of tau_l and tau_h
    ## basic values of C_f and F_sig
    PRM.PC_f=PRM.GRD_nx*0.05
    PRM.PF_sig=PRM.GRD_nx*0.05
    PRM.PC_n = 2; #2 fishers!
    
    JUMP = logspace(log10(0.001),log10(.1),32);
    CATCH  = logspace(log10(0.01),log10(1),32);
    result = cell(size(JUMP,1),size(CATCH,1));

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
            res =[PRM.PS_p,PRM.PC_q],cons.measure["Ts"],cons.measure["f1"],cons.measure["fij"],cons.measure["H"];

            println("With full information")
            #With information
            PRM.PC_lambda=1;
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            FLAGS["measure_frac"]=true
            
#            println(cons.SN)
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            Tsinf=cons.measure["Ts"]
            F1inf=cons.measure["f1"] 

            result[i,j,:] =tuple(res..., (cons.measure["Ts"],cons.measure["f1"],cons.measure["fij"],cons.measure["H"])...)

        end
    end
    xs = zeros(Float64,size(JUMP,1),size(CATCH,1),2)
    TS = zeros(Float64,size(JUMP,1),size(CATCH,1),PRM.PC_n)
    TSinf = zeros(Float64,size(JUMP,1),size(CATCH,1),PRM.PC_n)
    F1 = zeros(Float64,size(JUMP,1),size(CATCH,1),PRM.PC_n)
    F1inf = zeros(Float64,size(JUMP,1),size(CATCH,1),PRM.PC_n)
    Catch = zeros(Float64,size(JUMP,1),size(CATCH,1),PRM.PC_n)
    Catchinf = zeros(Float64,size(JUMP,1),size(CATCH,1),PRM.PC_n)
    Fij = zeros(Float64,size(JUMP,1),size(CATCH,1),PRM.PC_n)
    Fijinf = zeros(Float64,size(JUMP,1),size(CATCH,1),PRM.PC_n)
    for i = 1:size(result,1)
        for j=1:size(result,2)
            xs[i,j,:]= result[i,j][1]
            TS[i,j,:] = result[i,j][2]
            F1[i,j,:] = result[i,j][3]
            Fij[i,j,:] = result[i,j][4]
            Catch[i,j,:] = result[i,j][5]
            TSinf[i,j,:] = result[i,j][6]
            F1inf[i,j,:] = result[i,j][7]
            Fijinf[i,j,:] = result[i,j][8]
            Catchinf[i,j,:] = result[i,j][9]
        end
    end
    npzwrite("./Data/Data_Fig3_xs.npy", xs)
    npzwrite("./Data/Data_Fig3_noinf.npy", TS)
    npzwrite("./Data/Data_Fig3_inf.npy", TSinf)
    npzwrite("./Data/Data_Fig3_f1_noinf.npy", F1)
    npzwrite("./Data/Data_Fig3_f1_inf.npy", F1inf)
    npzwrite("./Data/Data_Fig3_catch_noinf.npy", Catch)
    npzwrite("./Data/Data_Fig3_catch_inf.npy", Catchinf)
    npzwrite("./Data/Data_Fig3_fij_noinf.npy", Fij)
    npzwrite("./Data/Data_Fig3_fij_inf.npy", Fijinf)
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
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS,3);
            result[i,j] = [cons.measure["Hrate"][1] cons.measure["Hdist"][1] cons.measure["f2"][1]];
            println("$i $j")
        end
    end

    TS = zeros(Float64,size(RP,1),size(NS,1),3)
    xs = zeros(Float64,size(RP,1),size(NS,1),2)
    for i = 1:length(RP)
        for j=1:length(NS)
            TS[i,j,:] = result[i,j]
            xs[i,j,:] = [RP[i]  NS[j] ]
        end
    end
    
    npzwrite("./Data/Data_Fig4opt.npy", TS)
    npzwrite("./Data/Data_Fig4opt_xs.npy", xs)
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
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS,2);
            result[i,j] = [cons.measure["Hrate"][1] cons.measure["Hdist"][1] cons.measure["f2"][1]];
            println("$i $j")
        end
    end

    TS = zeros(Float64,size(RP,1),size(NS,1),3)
    xs = zeros(Float64,size(RP,1),size(NS,1),2)
    for i = 1:length(RP)
        for j=1:length(NS)
            TS[i,j,:] = result[i,j]
            xs[i,j,:] = [RP[i]  NS[j] ]
        end
    end
    
    npzwrite("./Data/Data_Fig4opt_cliq.npy", TS)
    npzwrite("./Data/Data_Fig4opt_cliq_xs.npy", xs)
end

######## Test effect of friendship, on Tau_s and CPUE, as F_n is varied

##################### 1 FISH N FISHERS (do 5, 10, 20 fishers)
######## Stupid Optimal Social Network test (just keep adding links)

######## Evolving social network - find the ESS 








######## END  #######


