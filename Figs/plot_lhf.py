from plot_benchcore import *

pop=1 #LHF IFQ vs TAC total population x quota
lam=0 #LHF IFQ vs TAC lambda x quota
landscape=0 #LHF IFQ vs TAC  tauh x taul
depl=1 #LHF Hunting for conditions of local depletion effect


if pop:
    print "LHF IFQ vs TAC total population x quota"
    
    filename="Data_LHFpopIFQ"
    set_constants(filename)
    data_IFQ=load_data(filename) 
    filename="Data_LHFpopTAC"
    data_TAC=load_data(filename) 
    
    data=data_IFQ

    X=np.log10(data['pop'])
    Y=np.log10(data['quota'])
    TS=data['tausr']


    variables= ('Hrate','Hdist','Hrate_std','Hflux_std')
    labels=(r'$H$ rate','Hdist', r'std($H$ rate)',r'std(flux)')

    for i in range(len(variables)):
        H=data[variables[i] ]
        #H=np.log10(H)
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        plt.xlabel(r"Population")
        plt.ylabel(r"Quota")
        ax.set_zlabel(labels[i] )
        ax.plot_wireframe(X,Y,H)
        ax.scatter(X.ravel(),Y.ravel(),(H[:,:].ravel()))
        plt.show()
        
        if 0:
            mycmap = plt.cm.get_cmap('seismic')
            #mycmap.set_under('w')
            plt.xscale("log")
            plt.yscale("log")
            plt.xlim([np.min(X[X!=0]),np.max(X)])
            plt.ylim([np.min(Y),np.max(Y)])
            plt.pcolor(X,Y,H[:,:], cmap=mycmap,vmin =min(H[:,:].ravel()), vmax=max(H[:,:].ravel()))
            plt.colorbar()
            plt.show()



if depl:
    print "LHF IFQ vs TAC quota x PC_rp (Depletion effect)"
    
    filename="Data_LHFdeplIFQ"
    set_constants(filename)
    data_IFQ=load_data(filename) 
    filename="Data_LHFdeplTAC"
    data_TAC=load_data(filename) 
    
    data=data_IFQ

    X=(data['PC_rp'])
    Y=np.log10(data['quota']/constants["PF_n"])
    TS=data['tausr']


    variables= ('Hrate','Hdist','Hrate_std','Hflux_std')
    labels=(r'$H$ rate ','Hdist', r'std($H$ rate)',r'std(flux)')

    for i in range(len(variables)):
        H=data[variables[i] ]
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        plt.xlabel(r"PC_rp")
        plt.ylabel(r"Quota")
        ax.set_zlabel(labels[i] )
        ax.scatter(X.ravel(),Y.ravel(),np.log10(H.ravel()),c=Y.ravel())
        ax.plot_wireframe(X,Y,np.log10(H),color='grey')
        plt.show()
        #plot(*[l for l in H],xs=X[1,:])
        

if lam:
    #IFQ vs TAC, lambdaxquota
    filename="Data_IFQ"
    set_constants(filename)
    data_IFQ=load_data(filename) 
    filename="Data_TAC"
    data_TAC=load_data(filename) 
    
    data=data_IFQ
    
    X=data['PC_lambda']
    Y=np.log10(data['quota'])
    TS=data['tausr']
    #H=data['H']


    #Catch versus theory
    H=(data['Hrate'])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"quota")
    ax.set_zlabel(r"$H$")
    Hmax=np.zeros((X.shape[0],X.shape[1]))
    Hth=np.zeros((X.shape[0],X.shape[1]))
    TSth=np.zeros((X.shape[0],X.shape[1]))
    GRD_mx=constants["GRD_nx"]*constants["GRD_dx"]/2.
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            Hmax[i,j]=constants["PC_q"]*constants["PS_n"]*(constants["PC_f"]+constants["PF_sig"])**2 *pi/GRD_mx**2
            Hth[i,j]=fsbtheo(X[i,j],*get_constants(PC_lambda=X[i,j])  )[-2]*constants["PC_q"] #
            #TSth[i,j]=tausr(*get_constants() )
    
    ax.scatter(X.ravel(),Y.ravel(),H.ravel())
    if 0:
        filename2="Data_IFQ_cliq"
        #set_constants(filename)
        data_IFQc=load_data(filename2) 
        Hc=np.log10(data_IFQc['Hrate'])
        Xc=data_IFQc['PC_ncliq']
        Yc=np.log10(data_IFQc['quota'])
        ax.plot_wireframe(1/Xc,Yc,Hc,color='r')
    ax.scatter(X.ravel(),Y.ravel(),Hth.ravel(),color='g')
    ax.plot_wireframe(X,Y,Hth,color='g')
    ax.plot_wireframe(X,Y,Hmax,color='pink')
    plt.show()



    variables= ('H','Hrate','Hstd','Hfluxstd')
    labels=(r'$H$ dist',r'$H$ rate', r'std(H)',r'std(flux)')

    for i in range(len(variables)):
        H=data[variables[i] ]
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        plt.xlabel(r"$\lambda$")
        plt.ylabel(r"Quota")
        ax.set_zlabel(labels[i] )
        ax.scatter(X.ravel(),Y.ravel(),np.log10(H[:,:].ravel()))
        plt.show()
        
        if 0:
            mycmap = plt.cm.get_cmap('seismic')
            #mycmap.set_under('w')
            plt.xscale("log")
            plt.yscale("log")
            plt.xlim([np.min(X[X!=0]),np.max(X)])
            plt.ylim([np.min(Y),np.max(Y)])
            plt.pcolor(X,Y,H[:,:], cmap=mycmap,vmin =min(H[:,:].ravel()), vmax=max(H[:,:].ravel()))
            plt.colorbar()
            plt.show()


if 0:
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
    



if landscape:
    print 'LHF IFQ vs TAC:   tauh x taul'
    #IFQ vs TAC
    filename="Data_LHFvs2_IFQ"
    set_constants(filename)
    data_IFQ=load_data(filename) 
    filename="Data_LHFvs2_TAC"
    data_TAC=load_data(filename) 

    X=data_IFQ['PS_p']
    Y=data_IFQ['PC_q']
    taus1=X.copy()
    for i in range(taus1.shape[0]):
        for j in range(taus1.shape[1]):
            args=get_constants(PS_p=data_IFQ['PS_p'][i,j],PC_q=data_IFQ['PC_q'][i,j],PC_rp=data_IFQ['PC_rp'][i,j] )
            taus1[i,j]=tausr(*args)

    X=1./X
    Y=constants["PF_n"]/Y
    #X=data_IFQ['quota']
    X,Y=np.log10(X/taus1),np.log10(Y/taus1)

    

    T1,T2=(data_IFQ["Hrate_std"]).ravel(),(data_TAC["Hrate_std"]).ravel()
    #T1,T2=data_IFQ["Hflux_std"].ravel(),data_TAC["Hflux_std"][:,:][:,:].ravel()
    #T1,T2=data_IFQ["Hrate"].ravel(),data_TAC["Hrate"][:,:][:,:].ravel()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\log_{10}(\tau_l/\tau_s^{1})$")
    plt.ylabel(r"$\log_{10}(\tau_h/\tau_s^{1})$")
    ax.set_zlabel("Inequality (blue=TAC)" )
    Z=T2/T1-1
    Z=np.log10(Z)
    #ax.scatter(X.ravel(),Y.ravel(),Z,c=Z)
    ax.scatter(X.ravel(),Y.ravel(),np.log10(T1),c='r')
    ax.scatter(X.ravel(),Y.ravel(),np.log10(T2),c='b')
    
    #ax.set_zlim(bottom=min(np.minimum(T1,T2)), top=max(np.maximum(T1,T2)))
    #ax.set_zlim(bottom=min(Z), top=max(Z))
    #ax.set_zlabel(r"$stddev(H^{TAC})/stddev(H^{ITQ})-1$")
    plt.show()


