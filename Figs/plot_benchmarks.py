from plot_benchcore import *

#Plotter: numerical benchmarks
#see plot_bench_analytical for results that
#are significant only for the analytics


#### SWITCHES FOR WHICH FIGURES TO PLOT ####

firstpass=0
firstpass_ns=0
fig2a=0
fig2b=0
fig3=0 #Fig3 for tausr
fig3H=0 #Fig3 for catch rate
fig3f=0 #Other quantities related to fig3 (analytics)
fig4opt=0 #optimal lambda
fig4opt_cliq=0 #optimal ncliques
fig4opt_comp=0 #cliques vs lambda
rndcliq=0  #random partition of fishers into cliques
rndcliq_explor=1  #random partition of fishers into cliques for all tauh taul
spying=0
worst=1 # Look for worst value of lambda depending on tauh, taul


#=========== FIGURES ==================

if firstpass:
    filename="Data_firstpass"
    set_constants(filename)
    data=load_data(filename) 
    TS=data['\\tau_s^R']
    xs=data['PC_rp']
    #print(TS)
    plt.xlabel(r'$C_p$')
    plt.ylabel('First passage time')
    plot(TS,[tausr(*get_constants(PC_rp=x)) for x in xs],xs=xs)


if firstpass_ns:
    filename="Data_firstpass_ns"
    set_constants(filename)
    data=load_data(filename) 
    TS=data['\\tau_s^R']
    dist=data['dist']
    xs=data['PS_n']
    plt.xlabel(r'$S_n$')
    plt.ylabel('First passage time')
    plot(TS,[tausr(*get_constants(PS_n=x)) for x in xs],xs=xs,log='xy')
   # plt.plot(xs,[TS[0]/x for x in xs])
   # plot(xs,dist,hold=1)
    #plot(xs,[domain_size(*get_constants(PS_n=x)) for x in xs],log='xy')

if fig2a:
    #=========== FIGURE 2A ==================
    filename="Data_Fig2a"
    set_constants(filename)
    data=load_data(filename) 
    TS=data['\\tau_s^R']
    F1=data["f1"]
    xs=data['PC_rp']
    plt.plot(xs,TS)
    plt.plot(xs,[tausr(*get_constants(PC_rp=x)) for x in xs])
    
    plt.show()
    plt.plot(xs,F1)
    plt.plot(xs,[f1r(*get_constants(PC_rp=x)) for x in xs])
    plt.show()

if spying:
    #=========== FIGURE 2A ==================
    filename="Data_spy"
    set_constants(filename)
    data=load_data(filename) 
    TS=data['tau_s']
    H1=data['H1']
    H2=data['H2']
    xs=data['PC_spy']
    #plt.plot(xs,TS)
    #plt.show()
    plt.plot(xs,H1)
    plt.plot(xs,H2)
    plt.xlabel('Spying radius')
    plt.ylabel('Catch rate')
    plt.legend(["Sheep","Wolves"])
    plt.show()

if fig2b:
    #=========== FIGURE 2B ==================

    filename="Data_Fig2b"
    set_constants(filename)
    data=load_data(filename) 
    TS=data['\\tau_s^R']
    F1=data["f1"]
    X=data['PF_sig']
    Y=data['PC_f']
    Cp=data['PC_rp']

    #TS = np.load("../Data/Data_Fig2b.npy")
    #xy = np.load("../Data/Data_Fig2b_xs.npy")
    GRD_mx = constants["GRD_mx"]

    #3dscatter
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel("F_sigma")
    plt.ylabel("C_f")
    #plt.zlabel("tau_s")
    ax.set_zlim(bottom=0, top=500)
    XX,YY=X.ravel(),Y.ravel()
    ax.scatter(XX/GRD_mx,YY/GRD_mx,TS.ravel(),c=TS.ravel())

    #Wireframe
    T=X.copy()
    Texp=X.copy()
    for i in range(T.shape[0]):
        for j in range(T.shape[1]):
            args=get_constants(PF_sig=X[i,j],PC_f=Y[i,j],PC_rp= Cp[i,j])
            args[0]=p_opt(*args) #optimal turn probability
            T[i,j]=tausr(*args)
            args[0]=data["PC_rp"][i,j] #optimal turn probability
            Texp[i,j]=tausr(*args)


    #ax.plot_wireframe(X/GRD_mx,Y/GRD_mx,T)
    ax.plot_wireframe(X/GRD_mx,Y/GRD_mx,Texp)
    plt.show()

    #Heatmap
    Z=TS
    fig = plt.figure()
    plt.xlabel("F_sigma")
    plt.ylabel("C_f")
    ax = fig.add_subplot(111)
    ax.pcolor(X/GRD_mx, Y/GRD_mx, Z,  vmin=abs(Z).min(), vmax=500)#abs(Z).max())
    plt.show()

if fig3:
    #=========== FIGURE 3 ==================

    filename="Data_Fig3_inf"
    set_constants(filename)
    data_inf=load_data(filename) 

    filename="Data_Fig3_noinf"
    data_noinf=load_data(filename) 
    TS=data_noinf['\\tau_s^R']
    TSinf=data_inf['\\tau_s^R']
    X=data_inf['PS_p']
    Y=data_inf['PC_q']


    #3dscatter
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\log_{10}(\tau_l=1./S_p)/\tau_s^{1}$")
    plt.ylabel(r"$\log_{10}(\tau_h= F_n/C_q)/\tau_s^{1}$")
    #plt.zlabel("tau_s")
    Z1,Z2=TS[:,:].ravel(),TSinf[:,:].ravel()
    X=1./X
    Y=constants["PF_n"]/Y
    taus1=X.copy()
    for i in range(taus1.shape[0]):
        for j in range(taus1.shape[1]):
            args=get_constants(PS_p=data_inf['PS_p'][i,j],PC_q=data_inf['PC_q'][i,j],PC_rp=data_inf['PC_rp'][i,j] )
            taus1[i,j]=tausr(*args)
    X,Y=np.log10(X/taus1),np.log10(Y/taus1)

#    ax.set_zlim(bottom=-1, top=1)
    ax.scatter(X.ravel(),Y.ravel(),Z1/Z2-1,c=Z1/Z2-1)
    ax.set_zlabel(r"$VOI= \tau_s^R/\tau_s^I-1$")
    plt.show()
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\log_{10}(\tau_l=1./S_p)/\tau_s^{1}$")
    plt.ylabel(r"$\log_{10}(\tau_h= F_n/C_q)/\tau_s^{1}$")
    ax.set_zlim(bottom=0, top=200)
    ax.scatter(X.ravel(),Y.ravel(),Z1,c='r')
    ax.scatter(X.ravel(),Y.ravel(),Z2,c='b')
    ax.set_zlabel(r"$\tau_s$")

    #Wireframe
    T=X.copy()
    for i in range(T.shape[0]):
        for j in range(T.shape[1]):
            args=get_constants(PS_p=data_inf['PS_p'][i,j],PC_q=data_inf['PC_q'][i,j],PC_rp=data_inf['PC_rp'][i,j] )
            #args[0]=p_opt(*args) #optimal turn probability
            T[i,j]=tausr(*args)

    #ax.plot_wireframe(X,Y,T)
    #plt.show()

    #Heatmap
    #Z1=TS
    #Z2=TSinf
    #fig = plt.figure()
    #plt.xlabel(r"$\log_{10}\tau_l=1./S_p$")
    #plt.ylabel(r"$\log_{10}\tau_h= F_n/C_q$")
    #ax = fig.add_subplot(111)
    #ax.pcolor(X, Y, Z1-Z2,  vmin=(Z1-Z2).min(), vmax=(Z1-Z2).max())#abs(Z).max())
    #plt.show()
    
    
if fig3H:
    #CATCH
    filename="Data_Fig3_inf"
    set_constants(filename)
    data_inf=load_data(filename) 
    filename="Data_Fig3_noinf"
    data_noinf=load_data(filename) 
    Catch=data_noinf["H"]
    Catchinf=data_inf["H"]
    X=data_inf['PS_p']
    Y=data_inf['PC_q']
    taus1=X.copy()
    for i in range(taus1.shape[0]):
        for j in range(taus1.shape[1]):
            args=get_constants(PS_p=data_inf['PS_p'][i,j],PC_q=data_inf['PC_q'][i,j],PC_rp=data_inf['PC_rp'][i,j] )
            taus1[i,j]=tausr(*args)

    X=1./X
    Y=constants["PF_n"]/Y
    X,Y=np.log10(X/taus1),np.log10(Y/taus1)

    T1,T2=Catch[:,:].ravel(),Catchinf[:,:].ravel()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\log_{10}(\tau_l=1./S_p)/\tau_s^{1}$")
    plt.ylabel(r"$\log_{10}(\tau_h= F_n/C_q)/\tau_s^{1}$")
    Z=T2/T1-1
    ax.scatter(X.ravel(),Y.ravel(),Z,c=Z)
    ax.set_zlim(bottom=min(Z), top=max(Z))
    ax.set_zlabel(r"$H^I/H^R-1$")
    plt.show()
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\log_{10}\tau_l=1./S_p$")
    plt.ylabel(r"$\log_{10}\tau_h= F_n/C_q$")
    ax.scatter(X.ravel(),Y.ravel(),T1,c='r')
    ax.scatter(X.ravel(),Y.ravel(),T2,c='b')
    ax.set_zlabel(r"$H$")
    plt.show()
    
    mycmap = plt.cm.get_cmap('seismic')
    plt.xlabel(r"$\log_{10}\tau_l/\tau_s^{1}$")
    plt.ylabel(r"$\log_{10}\tau_h/\tau_s^{1}$")
    Z=Catchinf/Catch-1.
    plt.pcolor(X, Y, Z,  vmin=Z.min(), vmax=Z.max(),cmap=mycmap)
    plt.colorbar()
    plt.show()
    
if fig3f:
    #FRACTIONS
    filename="Data_Fig3_inf"
    set_constants(filename)
    data_inf=load_data(filename) 
    filename="Data_Fig3_noinf"
    data_noinf=load_data(filename) 
    TS=data_noinf['\\tau_s^R']
    TSinf=data_inf['\\tau_s^R']
    F1=data_noinf["f1"]
    F1inf=data_inf["f1"]
    Fij=data_noinf["fij"]
    Fijinf=data_inf["fij"]
    Catch=data_noinf["H"]

    X=data_inf['PS_p']
    Y=data_inf['PC_q']

    X=1./X
    Y=constants["PF_n"]/Y
    taus1=X.copy()
    for i in range(taus1.shape[0]):
        for j in range(taus1.shape[1]):
            args=get_constants(PS_p=data_inf['PS_p'][i,j],PC_q=data_inf['PC_q'][i,j],PC_rp=data_inf['PC_rp'][i,j] )
            taus1[i,j]=tausr(*args)

    X,Y=np.log10(X/taus1),np.log10(Y/taus1)
        
    #F1
    T1,T2=F1[:,:].ravel(),F1inf[:,:].ravel()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\log_{10}\tau_l=1./S_p$")
    plt.ylabel(r"$\log_{10}\tau_h= F_n/C_q$")
    ax.set_zlim(bottom=0, top=1)
    ax.scatter(X.ravel(),Y.ravel(),T1,c='r')
    ax.scatter(X.ravel(),Y.ravel(),T2,c='b')
    ax.set_zlabel(r"$f_1$")
    plt.show()
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\log_{10}\tau_l=1./S_p$")
    plt.ylabel(r"$\log_{10}\tau_h= F_n/C_q$")
    ax.set_zlim(bottom=-.5, top=.5)
    ax.scatter(X.ravel(),Y.ravel(),T1-T2,c=T1-T2)
    ax.set_zlabel(r"$f_1^R - f_1^I$")
    plt.show()
    
    #FIJ
    T1,T2=Fij[:,:].ravel(),Fijinf[:,:].ravel()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\log_{10}\tau_l=1./S_p$")
    plt.ylabel(r"$\log_{10}\tau_h= F_n/C_q$")
    ax.scatter(X.ravel(),Y.ravel(),T1-T2,c=T1-T2)
    ax.set_zlim(bottom=min(T1-T2), top=max(T1-T2))
    ax.set_zlabel(r"$f_{ij}^R - f_{ij}^I$")
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\log_{10}\tau_l=1./S_p$")
    plt.ylabel(r"$\log_{10}\tau_h= F_n/C_q$")
    ax.set_zlim(bottom=0, top=1)
    ax.scatter(X.ravel(),Y.ravel(),T1,c='r')
    ax.scatter(X.ravel(),Y.ravel(),T2,c='b')
    ax.set_zlabel(r"$f_{ij}$")
    plt.show()
    

if fig4opt:
    #=========== FIGURE 4 - OPTIMIZATION ==================
    filename="Data_Fig4opt"
    set_constants(filename)
    data=load_data(filename) 
    X=data['PC_lambda']
    Y=data['PS_n']
    TS=data['\\tau_s^R']
    F1=data['f1']
    F2=data['f2']
    H=data['Hrate']

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$S_n$")
#    ax.set_zlim(bottom=0, top=.01)
    ax.scatter(X.ravel(),Y.ravel(),TS[:,:].ravel())
    ax.set_zlabel(r"$\tau_s$")
    plt.show()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$S_n$")
    ax.set_zlabel(r"$H$")
    ax.scatter(X.ravel(),Y.ravel(),H[:,:].ravel())
    plt.show()
    
    mycmap = plt.cm.get_cmap('seismic')
    #mycmap.set_under('w')
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim([np.min(X[X!=0]),np.max(X)])
    plt.ylim([np.min(Y),np.max(Y)])
    plt.pcolor(X,Y,TS[:,:], cmap=mycmap,vmin =min(TS[:,:].ravel()), vmax=max(TS[:,:].ravel()))
    plt.colorbar()
    plt.show()
    

    
if fig4opt_cliq:
    filename="Data_Fig4opt_cliq"
    set_constants(filename)
    data=load_data(filename) 
    X=data['PC_ncliq']
    Y=data['PS_n']
    TS=data['\\tau_s^R']
    H=data['Hrate']

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$N_{cliques}$")
    plt.ylabel(r"$S_n$")
#    ax.set_zlim(bottom=0, top=.01)
    ax.scatter(X.ravel(),Y.ravel(),TS[:,:].ravel())
    ax.set_zlabel(r"$\tau_s^r$")
    plt.show()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$N_{cliques}$")
    plt.ylabel(r"$S_n$")
    ax.scatter(X.ravel(),Y.ravel(),H[:,:].ravel())
    ax.set_zlabel(r"$H$")
    plt.show()

    mycmap = plt.cm.get_cmap('seismic')
    #mycmap.set_under('w')
    #plt.xscale("log")
    #plt.yscale("log")
    plt.xlabel(r"$N_{cliques}$")
    plt.ylabel(r"$S_n$")
    plt.xlim([np.min(X),np.max(X)])
    plt.ylim([np.min(Y),np.max(Y)])
    plt.pcolor(X,Y,TS[:,:], cmap=mycmap,vmin =min(TS[:,:].ravel()), vmax=max(TS[:,:].ravel()))
    plt.colorbar()
    plt.show()
    
if fig4opt_comp:
    filename="Data_Fig4opt_comp"
    set_constants(filename)
    data=load_data(filename) 
    X=data['PC_ncliq']
    Y=data['PC_lambda']
    TS=data['\\tau_s^R']
    H=data['Hrate']

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$N_{cliques}$")
    plt.ylabel(r"$\lambda$")
#    ax.set_zlim(bottom=0, top=.01)
    ax.scatter(X.ravel(),Y.ravel(),TS[:,:].ravel())
    ax.set_zlabel(r"$\tau_s^r$")
    plt.show()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$N_{cliques}$")
    plt.ylabel(r"$\lambda$")
    ax.scatter(X.ravel(),Y.ravel(),H[:,:].ravel())
    ax.set_zlabel(r"$H$")
    plt.show()

    mycmap = plt.cm.get_cmap('Greys')
    #mycmap.set_under('w')
    #plt.xscale("log")
    #plt.yscale("log")
    plt.xlim([np.min(X),np.max(X)])
    plt.ylim([np.min(Y),np.max(Y)])
    plt.pcolor(X,Y,TS[:,:], cmap=mycmap,vmin =min(TS[:,:].ravel()), vmax=max(TS[:,:].ravel()))
    plt.colorbar()
    plt.show()
    
if rndcliq:
    filename="Data_rndcliq"
    set_constants(filename)
    data=load_data(filename) 
    Htyp="Hrate"
    occur=data["occur"]
    sizes=np.array([i for i,j in enumerate(occur) if j>1]) #clique sizes with more than one occurrence
    H=data[Htyp][sizes]
    Hvar=data[Htyp+"_var"][sizes]*occur[sizes]/(occur[sizes]-1) #unbiased variance
    Hstd=np.sqrt(Hvar)
    TSvar=data["TS_var"][sizes]*occur[sizes]/(occur[sizes]-1) #unbiased variance
    TSstd=np.sqrt(TSvar)
    TS=data["TS"][sizes]
    errorbar(sizes+1,H,yerr=Hstd)
    errorbar(sizes+1,TS,yerr=TSstd)
    
if rndcliq_explor:
    filename="Data_rndcliq_explor"
    set_constants(filename)
    data=load_data(filename) 
    Htyp="Hrate"
    occur=data["occur"]
    
    X=data['PS_p'][:,:,0]
    Y=data['PC_q'][:,:,0]
    Cp=data['PC_rp'][:,:,0]
    
    Z=np.zeros(occur.shape[:2])
    taus1=np.zeros(occur.shape[:2])
    for x in range(occur.shape[0]):
        for y in range(occur.shape[1]):
            loccur=occur[x,y]
            sizes=np.array([i for i,j in enumerate(loccur) if j>1]) #clique sizes with more than one occurrence
            H=data[Htyp][x,y,sizes]
            Hvar=data[Htyp+"_var"][x,y,sizes]*loccur[sizes]/(loccur[sizes]-1) #unbiased variance
            Hstd=np.sqrt(Hvar)
            TSvar=data["TS_var"][x,y,sizes]*loccur[sizes]/(loccur[sizes]-1) #unbiased variance
            TSstd=np.sqrt(TSvar)
            TS=data["TS"][x,y,sizes]
            Z[x,y]=1+sizes[np.argmax(H)]
            taus1[x,y]=tausr(*get_constants(PC_rp=Cp[x,y], PS_p=X[x,y], PC_q=Y[x,y] ))
            #errorbar(sizes+1,H,yerr=Hstd,hold=1)
#    errorbar(sizes+1,TS,yerr=TSstd)

    X=1./X
    Y=constants["PF_n"]/Y
    X,Y=np.log10(X/taus1),np.log10(Y/taus1)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\log_{10}\tau_l/\tau_s^{1}$")
    plt.ylabel(r"$\log_{10}\tau_h/\tau_s^{1}$")
    ax.scatter(X.ravel(),Y.ravel(),Z.ravel(),c=Z.ravel() )
    ax.set_zlabel(r"Optimal clique size")
    plt.show()
    plt.clf()
    mycmap = plt.cm.get_cmap('seismic')
    plt.xlabel(r"$\log_{10}\tau_l/\tau_s^{1}$")
    plt.ylabel(r"$\log_{10}\tau_h/\tau_s^{1}$")
    plt.pcolor(X, Y, Z,  vmin=Z.min(), vmax=Z.max(),cmap=mycmap)
    plt.colorbar()
    plt.show()

    

if worst:
    #CATCH
    filename="Data_worst"
    set_constants(filename)
    data=load_data(filename) 
    Catch=data["Hrate"]
    X=data['PS_p']
    Y=data['PC_q']
    L=data['PC_lambda']

    X=1./X
    Y=constants["PF_n"]/Y
    taus1=np.array([tausr(*get_constants(PC_rp=Cp )) for Cp in data['PC_rp']])
    X,Y=np.log10(X/taus1),np.log10(Y/taus1)

    XX=np.zeros((X.shape[0],X.shape[1]))
    YY=np.zeros((X.shape[0],X.shape[1]))
    Worst=np.zeros((X.shape[0],X.shape[1]))
    Dif=np.zeros((X.shape[0],X.shape[1]))
    Extrem=np.zeros((X.shape[0],X.shape[1]))
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            XX[i,j]=np.mean(X[i,j])
            YY[i,j]=np.mean(Y[i,j])
            Worst[i,j]=L[i,j,[z for z in range(X.shape[2]) if Catch[i,j,z]==min(Catch[i,j])][0]]
            Extrem[i,j]=(Catch[i,j,0]-Catch[i,j,-1])/min(Catch[i,j,0],Catch[i,j,1])
            Dif[i,j]=(max(Catch[i,j])/min(Catch[i,j]))-1.
            
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\log_{10}\tau_l/\tau_s^{1}$")
    plt.ylabel(r"$\log_{10}\tau_h/\tau_s^{1}$")
    ax.scatter(XX.ravel(),YY.ravel(),Worst.ravel(),c=Worst.ravel() )
    ax.set_zlabel(r"$Worst value of lambda$")
    plt.show()
    
    mycmap = plt.cm.get_cmap('seismic')
    plt.xlabel(r"$\log_{10}\tau_l/\tau_s^{1}$")
    plt.ylabel(r"$\log_{10}\tau_h/\tau_s^{1}$")
    Z=Worst
    plt.pcolor(XX, YY, Z,  vmin=Z.min(), vmax=Z.max(),cmap=mycmap)
    plt.colorbar()
    plt.show()

    print("H(best lambda)/H(worst lambda)-1")
    plt.clf()
    mycmap = plt.cm.get_cmap('seismic')
    plt.xlabel(r"$\log_{10}\tau_l/\tau_s^{1}$")
    plt.ylabel(r"$\log_{10}\tau_h/\tau_s^{1}$")
    Z=Dif
    plt.pcolor(XX, YY, Z,  vmin=Z.min(), vmax=Z.max(),cmap=mycmap)
    plt.colorbar()
    plt.show()
    
    print("H(best lambda)/H(worst lambda)-H(best lambda in 0 or 1)/H(worst lambda in 0 or 1) ")
    plt.clf()
    mycmap = plt.cm.get_cmap('seismic')
    plt.xlabel(r"$\log_{10}\tau_l/\tau_s^{1}$")
    plt.ylabel(r"$\log_{10}\tau_h/\tau_s^{1}$")
    Z=Dif-Extrem
    plt.pcolor(XX, YY, Z,  vmin=Z.min(), vmax=Z.max(),cmap=mycmap)
    plt.colorbar()
    plt.show()

