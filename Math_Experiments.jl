
##########################!!!!! SIMULATION EXPERIMENTS !!!!!!!!################

#Experiments geared toward Matthieu's math


function ana_flux()
    println("Measure fluxes between states in ABM")
    PRM.PC_f=PRM.GRD_mx*0.01
    PRM.PF_sig=PRM.GRD_mx*0.04
    PRM.PC_n = 10 #10 fishers
    PRM.PF_n = 100
    PRM.PC_lambda= 1.
    PRM.PS_n = int(round(.05 / (pi *PRM.PF_sig^2 /PRM.GRD_nx^2 ) )) ;
    
    #Tmptest: constants from fig4opt
    PRM.PC_ncliq = 1;
    PRM.PS_p = 0.1;
    PRM.PC_q = 0.1;
    PRM.PC_rp = fnc_optimal_PCrp();
    #/Tmptest

    
    JUMP = logspace(log10(0.001),log10(.01),4);
    CATCH  = logspace(log10(0.1),log10(1),4);
    LAM=linspace(0.01,1.,4);
    #FN=[2,4,8,12];
    #RP = logspace(log10(PRM.PF_n*PRM.PC_n/2),log10(PRM.PF_n*20*PRM.PC_n),4);
    RM=  logspace(log10(0.005),log10(0.02),4)*PRM.GRD_mx;
    result = cell(size(JUMP,1),size(CATCH,1));
    
    nens =40 #number of iterations for ensemble average
    
    headers=R"PS_p PC_q quota PC_rp PC_lambda tausr Hdist Hdist_std Hflux_std Htot Hrate Hrate_std"

    for i = 1:length(JUMP)
        for j = 1:length(CATCH)
            PRM.PC_f = RM[i]
            #PRM.PS_p = JUMP[i];
            #PRM.PC_q = CATCH[j];
            PRM.PC_lambda=LAM[j];
            #PRM.PF_n=FN[j]
            PRM.PC_rp = fnc_optimal_PCrp();
            
            tauh=PRM.PF_n/PRM.PC_q
            taul=1./PRM.PS_p
            taus=fnc_taus1()
            print("\n$i $j tau_h/tau_s: $(tauh/taus) tau_l/tau_s:  $(taul/taus)")
             
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            init_network(cons,FLAGS)
            quota=FLAGS["stopamount"]=convert(Int,PF_n*PS_n.*PC_n*200)
            FLAGS["measure_frac"]=true
            FLAGS["measure_flux"]=true
            FLAGS["stop"]=1
            #make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            cons.measure,cons.series=make_ensemble(nens,FLAGS)
            cons.measure["turn"]*=nens
            result[i,j] = JUMP[i], CATCH[j],quota,PRM.PC_rp,PRM.PC_lambda,mean(cons.measure["Ts"]),mean(cons.measure["Hdist"]),std(cons.measure["Hdist"]),std(cons.series["Hflux"]),mean(cons.measure["Htot"]),mean(cons.measure["Hrate"]),std(cons.measure["Hrate"]) ;
            for evt in ("found_school","left_school","lost_school", "bound","states","dist")
                npzwrite("./Data/$evt-$i-$j.npy", cons.series[evt] )
            end
            fout=open("./Data/turn$i-$j.dat", "w")
            turn=cons.measure["turn"]
            write(fout, "$(turn)\n" )
            taus=mean(cons.measure["Ts"])
            write(fout, "$(taus)\n" )
            taus=fnc_taus1()
            write(fout, "$(taus)\n" )
            for f in names(PRM)
                write(fout,"$f    $(getfield(PRM,f))\n" )
            end
            close(fout)
        end
    end
    
    save_results(result,headers,"Data_anaflux",PRM)

end
