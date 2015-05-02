
##########################!!!!! SIMULATION EXPERIMENTS !!!!!!!!################


function do_LHFmovie()
    IFQ=true #IFQ or TAC?

    PRM.PC_f=PRM.GRD_mx*0.03
    PRM.PF_sig=PRM.GRD_mx*0.04
    PRM.PC_n = 10 #10 fishers
    PRM.PF_n = 100
    PRM.PS_n = int(round(.05 / (pi *PRM.PF_sig^2 /PRM.GRD_nx^2 ) )) ;
    
    RP = logspace(log10(PRM.PF_n/2),log10(PRM.PF_n*PRM.PS_n.*10),8);

    PRM.PC_rp = fnc_optimal_PCrp();
    quota=RP[3]
    PRM.PC_lambda=0.

    school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
    init_network(cons,FLAGS)
    FLAGS["save"]=true
    FLAGS["savepath"]="./LHFMovie2/"
    FLAGS["measure_frac"]=false
    if IFQ
        FLAGS["stop"]=5
        FLAGS["IFQ"]=true
        FLAGS["stopamount"]=convert(Int,floor(quota))
    else:
        FLAGS["stopamount"]=convert(Int,floor(quota*PRM.PC_n))
        FLAGS["stop"]=1
    end
    make_season(school,fish,cons,fishtree,EVENTS,FLAGS,OUT);
    println("""Hrate $(mean(cons.measure["Hrate"]))  Hdist $(mean(cons.measure["Hdist"]))""")

end


## Clique size vs Lambda => equity, flux?


function do_TAC()
    ######  Vary the total quota
    println("TAC")
    PRM.PC_f=PRM.GRD_mx*0.03
    PRM.PF_sig=PRM.GRD_mx*0.04
    PRM.PC_n = 10 #10 fishers
    PRM.PF_n = 100
    PRM.PS_n = int(round(.05 / (pi *PRM.PF_sig^2 /PRM.GRD_nx^2 ) )) ;
    
    RP = logspace(log10(PRM.PF_n/2),log10(PRM.PF_n*PRM.PS_n.*10),8);
    LAM = linspace(0.,1.,8);
    result = cell(size(RP,1),size(LAM,1));
    PRM.PC_rp = fnc_optimal_PCrp();
    
    for i = 1:length(RP)
        for j=1:length(LAM)
            PRM.PC_lambda = LAM[j];
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            init_network(cons,FLAGS)
            FLAGS["measure_frac"]=false
            FLAGS["stopamount"]=convert(Int,floor(RP[i]*PRM.PC_n))
            FLAGS["stop"]=1
            #make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            nens=convert(Int,maximum([1,round(4000/log10(RP[i])^2)]) )
            println("$i $j $nens")
            cons.measure,cons.series=make_ensemble(nens,FLAGS)
            result[i,j] = RP[i], LAM[j],mean(cons.measure["Ts"]),mean(cons.measure["Hdist"]),std(cons.measure["Hdist"]),std(cons.series["Hflux"]),mean(cons.measure["Htot"]),mean(cons.measure["Hrate"]) ;
        end
    end
    save_results(result,R"quota PC_lambda tausr Hdist Hdist_std Hfluxstd Htot Hrate","Data_TAC",PRM)
end

function do_IFQ()
    ######  Vary the individual quota
    println("IFQ")
    PRM.PC_f=PRM.GRD_mx*0.03
    PRM.PF_sig=PRM.GRD_mx*0.04
    PRM.PC_n = 10 #10 fishers
    PRM.PF_n = 100
    PRM.PS_n = int(round(.05 / (pi *PRM.PF_sig^2 /PRM.GRD_nx^2 ) )) ;
    
    RP = logspace(log10(PRM.PF_n/2),log10(PRM.PF_n*PRM.PS_n.*10),8);
    LAM = linspace(0.,1.,8);
    result = cell(size(RP,1),size(LAM,1));
    PRM.PC_rp = fnc_optimal_PCrp();
    
    for i = 1:length(RP)
        for j=1:length(LAM)
            PRM.PC_lambda = LAM[j];
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            init_network(cons,FLAGS)
            FLAGS["measure_frac"]=false
            FLAGS["stopamount"]= convert(Int,round(RP[i]))
            FLAGS["stop"]=5
            FLAGS["IFQ"]=true
            #make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            nens=convert(Int,maximum([1,round(4000/log10(RP[i])^2)]) )
            println("$i $j $nens")
            cons.measure,cons.series=make_ensemble(nens,FLAGS)
            result[i,j] =FLAGS["stopamount"], LAM[j],mean(cons.measure["Ts"]),mean(cons.measure["Hdist"]),std(cons.measure["Hdist"]),std(cons.series["Hflux"]),mean(cons.measure["Htot"]),mean(cons.measure["Hrate"]) ;
        end
    end
    save_results(result,R"quota PC_lambda tausr Hdist Hdist_std Hfluxstd Htot Hrate","Data_IFQ",PRM)
end


function do_LHFpop()
    println("LHF: Population")
    IFQ=true
    PRM.PC_f=PRM.GRD_mx*0.03
    PRM.PF_sig=PRM.GRD_mx*0.04
    PRM.PC_n = 10 #10 fishers
    PRM.PF_n = 100
    PRM.PS_n = int(round(.05 / (pi *PRM.PF_sig^2 /PRM.GRD_nx^2 ) )) ;
    PRM.PS_p = 0.001
    
    RP = logspace(log10(5),log10(PRM.PF_n*PRM.PS_n./10),8);
    POP = logspace(log10(.1),log10(10),8)*PRM.PF_n*PRM.PS_n;
    PRM.PC_lambda = 1.;
    result = cell(size(RP,1),size(POP,1));
    
    for i = 1:length(RP)
        for j=1:length(POP)
            quota=iround(RP[i]*PRM.PC_n)
            pop=POP[j]
            #PRM.PF_n=maximum([1,iround(pop/PRM.PS_n)])
            PRM.PS_n=maximum([1,iround(pop/PRM.PF_n)])
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();

            PRM.PC_rp = fnc_optimal_PCrp();    
            init_network(cons,FLAGS)
            FLAGS["measure_frac"]=false
            if IFQ
                FLAGS["stopamount"]=iround(quota/PRM.PC_n)
                FLAGS["stop"]=5
                FLAGS["IFQ"]=true
             else
                FLAGS["stopamount"]=quota
                FLAGS["stop"]=1
            end
            #make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            nens=convert(Int,maximum([1,round(10/log10(RP[i])^2)]) )
            println("$i $j $nens")
            cons.measure,cons.series=make_ensemble(nens,FLAGS)
            result[i,j] = quota,pop,mean(cons.measure["Ts"]),mean(cons.measure["Hdist"]),cons.measure["Hdist_std"],std(cons.series["Hflux"]),mean(cons.measure["Hrate"]),cons.measure["Hrate_std"] ;
        end 
    end
    
    filename="Data_LHFpopSN"
    if IFQ
        filename="$(filename)IFQ"
    else
        filename="$(filename)TAC"
    end
    save_results(result,R"quota pop tausr Hdist Hdist_std Hflux_std Hrate Hrate_std",filename,PRM)

end


function do_LHFdepletion()
    println("LHF: Depletion effect vs PC_rp or PS_n")
    IFQ=true
    PRM.PC_f=PRM.GRD_mx*0.03
    PRM.PF_sig=PRM.GRD_mx*0.04
    PRM.PC_n = 10 #10 fishers
    PRM.PF_n = 100
    PRM.PS_n = int(round(.05 / (pi *PRM.PF_sig^2 /PRM.GRD_nx^2 ) )) ;
    PRM.PS_p = 0.001
    
    QUOTA = logspace(log10(5),log10(PRM.PF_n*PRM.PS_n./10),8);
    RP= logspace(log10(.25),log10(1.2),8)*fnc_optimal_PCrp();
    SP= logspace(log10(.1),log10(200),8)*0.001;
    PRM.PC_lambda = 1.;
    result = cell(size(QUOTA,1),size(RP,1));
    
    for i = 1:length(QUOTA)
        for j=1:length(SP)
            quota=iround(QUOTA[i]*PRM.PC_n)
            PRM.PC_rp=RP[j]
            #PRM.PS_p = SP[j]
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            init_network(cons,FLAGS)
            FLAGS["measure_frac"]=false
            if IFQ
                FLAGS["stopamount"]=iround(quota/PRM.PC_n)
                FLAGS["stop"]=5
                FLAGS["IFQ"]=true
             else
                FLAGS["stopamount"]=quota
                FLAGS["stop"]=1
            end
            nens=convert(Int,maximum([1,round(5000/log10(QUOTA[i])^2)]) )
            println("$i $j $nens")
            cons.measure,cons.series=make_ensemble(nens,FLAGS)
            result[i,j] = quota,PRM.PC_rp,PRM.PS_p,mean(cons.measure["Ts"]),mean(cons.measure["Hdist"]),cons.measure["Hdist_std"],std(cons.series["Hflux"]),mean(cons.measure["Hrate"]),cons.measure["Hrate_std"] ;
        end 
    end
    
    filename="Data_LHFdepl"
    if IFQ
        filename="$(filename)IFQ"
    else
        filename="$(filename)TAC"
    end
    save_results(result,R"quota PC_rp PS_p tausr Hdist Hdist_std Hflux_std Hrate Hrate_std",filename,PRM)

end


function do_IFQ_cliq()
    ######  Vary the individual quota
    println("IFQ_cliq")
    PRM.PC_f=PRM.GRD_mx*0.03
    PRM.PF_sig=PRM.GRD_mx*0.04
    PRM.PC_n = 10 #10 fishers
    PRM.PF_n = 100
    PRM.PC_lambda=1.
    PRM.PS_n = int(round(.05 / (pi *PRM.PF_sig^2 /PRM.GRD_nx^2 ) )) ;
    
    RP = logspace(log10(PRM.PF_n/2),log10(PRM.PF_n*PRM.PS_n.*10),8);
    LAM = [1,2,3,5,10];
    result = cell(size(RP,1),size(LAM,1));
    PRM.PC_rp = fnc_optimal_PCrp();
    
    for i = 1:length(RP)
        for j=1:length(LAM)
            PRM.PC_ncliq = convert(Int,round(LAM[j]));
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            init_network(cons,FLAGS)
            FLAGS["measure_frac"]=false
            FLAGS["stopamount"]= convert(Int,round(RP[i]))
            FLAGS["stop"]=5
            FLAGS["IFQ"]=true
            #make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            nens=convert(Int,maximum([1,round(500/log10(RP[i])^2)]) )
            println("$i $j $nens")
            cons.measure,cons.series=make_ensemble(nens,FLAGS)
            result[i,j] =FLAGS["stopamount"], PRM.PC_ncliq ,mean(cons.measure["Ts"]),mean(cons.measure["Hdist"]),std(cons.measure["Hdist"]),std(cons.series["Hflux"]),mean(cons.measure["Htot"]),mean(cons.measure["Hrate"]) ;
        end
    end
    save_results(result,R"quota PC_ncliq tausr Hdist Hdist_std Hfluxstd Htot Hrate","Data_IFQ_cliq",PRM)
end



function do_LHFvs()
    println("IFQ vs TAC ")
    PRM.PC_f=PRM.GRD_mx*0.03
    PRM.PF_sig=PRM.GRD_mx*0.04
    PRM.PC_n = 10 #10 fishers
    PRM.PF_n = 100
    PRM.PC_lambda= 1.
    PRM.PS_n = int(round(.05 / (pi *PRM.PF_sig^2 /PRM.GRD_nx^2 ) )) ;
    
    JUMP = logspace(log10(0.001),log10(.1),4);
    CATCH  = logspace(log10(0.1),log10(1),4);
    LAM=linspace(0.01,1.,4);
    RP = logspace(log10(PRM.PF_n*PRM.PC_n/2),log10(PRM.PF_n*20*PRM.PC_n),4);
    result_IFQ = cell(size(JUMP,1),size(CATCH,1));
    result_TAC = cell(size(JUMP,1),size(CATCH,1));
    
    nens =50 #number of iterations for ensemble average

    for i = 1:length(JUMP)
        for j = 1:length(CATCH)
            #PRM.PS_p = JUMP[i];
            #PRM.PC_lambda= LAM[i];
            PRM.PC_q = CATCH[j];
            quota=  round(RP[i])#PF_n*PC_n*5
            PRM.PC_rp = fnc_optimal_PCrp();
            print("\n$i $j ")
            
            print("TAC ")
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            init_network(cons,FLAGS)
            FLAGS["measure_frac"]=false
            FLAGS["stopamount"]= convert(Int,quota)
            FLAGS["stop"]=1
            #make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            cons.measure,cons.series=make_ensemble(nens,FLAGS)
            result_TAC[i,j] = JUMP[i], CATCH[j],quota,PRM.PC_rp,PRM.PC_lambda,mean(cons.measure["Ts"]),mean(cons.measure["Hdist"]),std(cons.measure["Hdist"]),std(cons.series["Hflux"]),mean(cons.measure["Htot"]),mean(cons.measure["Hrate"]),std(cons.measure["Hrate"]) ;
            
            print("IFQ\n")
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            init_network(cons,FLAGS)
            FLAGS["measure_frac"]=false
            FLAGS["stopamount"]= convert(Int,round(quota/PC_n))
            FLAGS["stop"]=5
            FLAGS["IFQ"]=true
            #make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            cons.measure,cons.series=make_ensemble(nens,FLAGS)
            result_IFQ[i,j] = JUMP[i], CATCH[j],quota,PRM.PC_rp,PRM.PC_lambda,mean(cons.measure["Ts"]),mean(cons.measure["Hdist"]),std(cons.measure["Hdist"]),std(cons.series["Hflux"]),mean(cons.measure["Htot"]),mean(cons.measure["Hrate"]),std(cons.measure["Hrate"]) ;

        end
    end
    save_results(result_TAC,R"PS_p PC_q quota PC_rp PC_lambda tausr Hdist Hdist_std Hflux_std Htot Hrate Hrate_std","Data_LHFvs_TAC",PRM)
    save_results(result_IFQ,R"PS_p PC_q quota PC_rp PC_lambda tausr Hdist Hdist_std Hflux_std Htot Hrate Hrate_std","Data_LHFvs_IFQ",PRM)

end



function do_LHFvs2()
    println("IFQ vs TAC 2")
    PRM.PC_f=PRM.GRD_mx*0.03
    PRM.PF_sig=PRM.GRD_mx*0.04
    PRM.PC_n = 10 #10 fishers
    PRM.PF_n = 100
    PRM.PC_lambda= 1.
    PRM.PS_n = int(round(.05 / (pi *PRM.PF_sig^2 /PRM.GRD_nx^2 ) )) ;
    
    JUMP = logspace(log10(0.001),log10(.1),4);
    CATCH  = logspace(log10(0.1),log10(1),4);
    LAM=linspace(0.01,1.,4);
    result_IFQ = cell(size(JUMP,1),size(CATCH,1));
    result_TAC = cell(size(JUMP,1),size(CATCH,1));
    
    quota=PRM.PF_n*PRM.PC_n
    
    nens =100 #number of iterations for ensemble average

    for i = 1:length(JUMP)
        for j = 1:length(CATCH)
            PRM.PS_p = JUMP[i];
            #PRM.PC_lambda= LAM[i];
            PRM.PC_q = CATCH[j];
            PRM.PC_rp = fnc_optimal_PCrp();
            print("\n$i $j ")
            
            print("TAC ")
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            init_network(cons,FLAGS)
            FLAGS["measure_frac"]=false
            FLAGS["stopamount"]= convert(Int,quota)
            FLAGS["stop"]=1
            #make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            cons.measure,cons.series=make_ensemble(nens,FLAGS)
            result_TAC[i,j] = JUMP[i], CATCH[j],quota,PRM.PC_rp,PRM.PC_lambda,mean(cons.measure["Ts"]),mean(cons.measure["Hdist"]),std(cons.measure["Hdist"]),std(cons.series["Hflux"]),mean(cons.measure["Htot"]),mean(cons.measure["Hrate"]),std(cons.measure["Hrate"]) ;
            
            print("IFQ\n")
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            init_network(cons,FLAGS)
            FLAGS["measure_frac"]=false
            FLAGS["stopamount"]= iround(quota/PC_n)
            FLAGS["stop"]=5
            FLAGS["IFQ"]=true
            #make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            cons.measure,cons.series=make_ensemble(nens,FLAGS)
            result_IFQ[i,j] = JUMP[i], CATCH[j],quota,PRM.PC_rp,PRM.PC_lambda,mean(cons.measure["Ts"]),mean(cons.measure["Hdist"]),std(cons.measure["Hdist"]),std(cons.series["Hflux"]),mean(cons.measure["Htot"]),mean(cons.measure["Hrate"]),std(cons.measure["Hrate"]) ;

        end
    end
    save_results(result_TAC,R"PS_p PC_q quota PC_rp PC_lambda tausr Hdist Hdist_std Hflux_std Htot Hrate Hrate_std","Data_LHFvs2_TAC",PRM)
    save_results(result_IFQ,R"PS_p PC_q quota PC_rp PC_lambda tausr Hdist Hdist_std Hflux_std Htot Hrate Hrate_std","Data_LHFvs2_IFQ",PRM)

end
