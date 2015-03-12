
##########################!!!!! SIMULATION EXPERIMENTS !!!!!!!!################

#Experiments geared toward Matthieu's math




function ana_flux()
    println("Measure fluxes between states in ABM")
    PRM.PC_f=PRM.GRD_mx*0.03
    PRM.PF_sig=PRM.GRD_mx*0.04
    PRM.PC_n = 10 #10 fishers
    PRM.PF_n = 100
    PRM.PC_lambda= 1.
    PRM.PS_n = int(round(.05 / (pi *PRM.PF_sig^2 /PRM.GRD_nx^2 ) )) ;
    
    JUMP = logspace(log10(0.001),log10(.1),5);
    CATCH  = logspace(log10(0.1),log10(1),5);
    LAM=linspace(0.01,1.,4);
    RP = logspace(log10(PRM.PF_n*PRM.PC_n/2),log10(PRM.PF_n*20*PRM.PC_n),4);
    result_IFQ = cell(size(JUMP,1),size(CATCH,1));
    result_TAC = cell(size(JUMP,1),size(CATCH,1));
    
    nens =50 #number of iterations for ensemble average

    for i = 1:length(JUMP)
        for j = 1:length(CATCH)
            PRM.PS_p = JUMP[i];
            PRM.PC_q = CATCH[j];
            PRM.PC_rp = fnc_optimal_PCrp();
            quota=FLAGS["stopamount"]
            print("\n$i $j ")
            
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            init_network(cons,FLAGS)
            FLAGS["measure_frac"]=true
            FLAGS["measure_flux"]=true
            FLAGS["stop"]=1
            #make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            cons.measure,cons.series=make_ensemble(nens,FLAGS)
            result_TAC[i,j] = JUMP[i], CATCH[j],quota,PRM.PC_rp,PRM.PC_lambda,mean(cons.measure["Ts"]),mean(cons.measure["Hdist"]),std(cons.measure["Hdist"]),std(cons.series["Hflux"]),mean(cons.measure["Htot"]),mean(cons.measure["Hrate"]),std(cons.measure["Hrate"]) ;
            

        end
    end
    save_results(result,R"PS_p PC_q quota PC_rp PC_lambda tausr Hdist Hdist_std Hflux_std Htot Hrate Hrate_std","Data_LHFvs_TAC",PRM)

end
