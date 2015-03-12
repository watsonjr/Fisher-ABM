
##########################!!!!! SIMULATION EXPERIMENTS !!!!!!!!################


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
    
    for i = 1:length(RP)
        for j=1:length(LAM)
            PRM.PC_lambda = LAM[j];
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            init_network(cons,FLAGS)
            FLAGS["measure_frac"]=false
            FLAGS["stopamount"]=convert(Int,floor(RP[i]))
            FLAGS["stop"]=1
            #make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            nens=convert(Int,maximum([1,round(4000/log10(RP[i])^2)]) )
            println("$i $j $nens")
            cons.measure,cons.series=make_ensemble(nens,FLAGS)
            result[i,j] = RP[i], LAM[j],mean(cons.measure["Ts"]),mean(cons.measure["Hdist"]),std(cons.measure["Hdist"]),std(cons.series["Hflux"]),mean(cons.measure["Htot"]),mean(cons.measure["Hrate"]) ;
        end
    end
    save_results(result,R"quota PC_lambda tausr H Hstd Hfluxstd Htot Hrate","Data_IFQ",PRM)
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
            #quota=  round(RP[i])#PF_n*PC_n*5
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

