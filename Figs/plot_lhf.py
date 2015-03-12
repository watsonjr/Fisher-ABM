from plot_benchcore import *

IFQ=0
LHFvs=1

if IFQ:
    #=========== FIGURE 4 - OPTIMIZATION ==================
    filename="Data_IFQ"
    set_constants(filename)
    data=load_data(filename) 
    X=data['PC_lambda']
    Y=np.log10(data['quota'])
    TS=data['tausr']
    H=data['H']

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


if LHFvs:
    #IFQ vs TAC
    filename="Data_LHFvs_IFQ"
    set_constants(filename)
    data_IFQ=load_data(filename) 
    filename="Data_LHFvs_TAC"
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
    X=data_IFQ['quota']
    X,Y=np.log10(X/taus1),np.log10(Y/taus1)


    T1,T2=(data_IFQ["Hrate_std"]/data_IFQ["Hrate"]).ravel(),(data_TAC["Hrate_std"]/data_TAC["Hrate"]).ravel()
    #T1,T2=data_IFQ["Hflux_std"].ravel(),data_TAC["Hflux_std"][:,:][:,:].ravel()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\log_{10}(\tau_l=1./S_p)/\tau_s^{1}$")
    plt.ylabel(r"$\log_{10}(\tau_h= F_n/C_q)/\tau_s^{1}$")
    plt.xlabel(r"Quota")
    plt.ylabel(r"Catchability")
    Z=T2/T1-1
    Z=np.log10(Z)
    #ax.scatter(X.ravel(),Y.ravel(),Z,c=Z)
    ax.scatter(X.ravel(),Y.ravel(),np.log10(T1),c='r')
    ax.scatter(X.ravel(),Y.ravel(),np.log10(T2),c='b')
    
    #ax.set_zlim(bottom=min(np.minimum(T1,T2)), top=max(np.maximum(T1,T2)))
    #ax.set_zlim(bottom=min(Z), top=max(Z))
    #ax.set_zlabel(r"$stddev(H^{TAC})/stddev(H^{ITQ})-1$")
    plt.show()


