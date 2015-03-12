from plot_benchcore import *
import random_partition as part

fig2b=0
fig2c=0
fig3=0
fig4opt =1
fig4opt_cn =0


    
#=========================== FIGURES =================================


set_constants("Data_params")



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
    GRD_mx = constants["GRD_mx"]


    #f1 versus theory
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$F_\sigma$")
    plt.ylabel(r"$C_f$")
    ax.set_zlabel(r"$f1$")
    F1th=np.zeros((X.shape[0],X.shape[1]))
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            F1th[i,j]=f1r(*get_constants(PF_sig=X[i,j],PC_f=Y[i,j],PC_rp= Cp[i,j]))
    ax.scatter(X.ravel(),Y.ravel(),F1[:,:].ravel())
    ax.plot_wireframe(X,Y,F1th[:,:],alpha=.5,color='r')
    #ax.scatter(X.ravel(),Y.ravel(),F1th[:,:].ravel(),color='r')
    plt.show()
    
if fig2c:
    #=========== FIGURE 2C ==================
    filename="Data_Fig2c"
    set_constants(filename)
    data=load_data(filename) 
    TS=data['\\tau_s^R']
    F1=data["f1"]
    X=data['S_p']
    Y=data['C_q']
    Cp=data['PC_rp']

    GRD_mx = constants["GRD_mx"]

    #f1 versus theory
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$S_p$")
    plt.ylabel(r"$C_q$")
    ax.xaxis.set_scale('log')
    ax.yaxis.set_scale('log')
    ax.set_xlim(min(X.ravel()),max(X.ravel()) )
    ax.set_ylim(min(Y.ravel()),max(Y.ravel()) )
    ax.set_zlabel(r"$f1$")
    F1th=np.zeros((X.shape[0],X.shape[1]))
    
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            F1th[i,j]=f1r(*get_constants(PS_p=X[i,j],PC_q=Y[i,j],PC_rp= Cp[i,j] ))
    ax.scatter(X.ravel(),Y.ravel(),F1[:,:].ravel())
    ax.plot_wireframe(X,Y,F1th[:,:],alpha=.5,color='r')
    #ax.scatter(X.ravel(),Y.ravel(),F1th[:,:].ravel(),color='r')
    plt.show()

    fig = plt.figure()
    plt.xlabel(r"$S_p$")
    plt.ylabel(r"$C_q$")
    ax = fig.add_subplot(111)
    ax.xaxis.set_scale('log')
    ax.yaxis.set_scale('log')
    ax.set_xlim(min(X.ravel()),max(X.ravel()) )
    ax.set_ylim(min(Y.ravel()),max(Y.ravel()) )
    plt.pcolor(X,Y,F1/F1th-1 )
    plt.colorbar()
    plt.show()



if fig3:
    #=========== FIGURE 3 ==================

    filename="Data_Fig3_inf"
    set_constants(filename)
    data_inf=load_data(filename) 

    filename="Data_Fig3_noinf"
    data_noinf=load_data(filename) 
    TS=data_noinf['fij']
    TSinf=data_inf['fij']
    X=data_inf['PS_p']
    Y=data_inf['PC_q']


    #3dscatter
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\log_{10}(\tau_l=1./S_p)/\tau_s^{1}$")
    plt.ylabel(r"$\log_{10}(\tau_h= F_n/C_q)/\tau_s^{1}$")
    Z1,Z2=TS[:,:].ravel(),TSinf[:,:].ravel()
    X=1./X
    Y=constants["PF_n"]/Y
    taus1=X.copy()
    for i in range(taus1.shape[0]):
        for j in range(taus1.shape[1]):
            args=get_constants(PS_p=data_inf['PS_p'][i,j],PC_q=data_inf['PC_q'][i,j],PC_rp=data_inf['PC_rp'][i,j] )
            taus1[i,j]=tausr(*args)
    X,Y=np.log10(X/taus1),np.log10(Y/taus1)
    ax.scatter(X.ravel(),Y.ravel(),Z2)
    ax.set_zlabel(r"$f(d_{ij}<\xi)$")

    #Wireframe
    B=X.copy()
    F1=X.copy()
    F2=X.copy()
    T=X.copy()
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            args=get_constants(PS_p=data_inf['PS_p'][i,j],PC_q=data_inf['PC_q'][i,j],PC_rp=data_inf['PC_rp'][i,j] )
            #args[0]=p_opt(*args) #optimal turn probability
            B[i,j]=fsbtheo(1.,*args)[3]
            F1[i,j]=fsbtheo(1.,*args)[0]
            F2[i,j]=fsbtheo(1.,*args)[1]
            T[i,j]=tausr(*args)
    ax.plot_wireframe(X,Y,B )
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_zlabel(r"$f_1$")
    plt.xlabel(r"$\log_{10}(\tau_l=1./S_p)/\tau_s^{1}$")
    plt.ylabel(r"$\log_{10}(\tau_h= F_n/C_q)/\tau_s^{1}$")
    ax.scatter(X.ravel(),Y.ravel(),data_inf["f1"] )
    ax.plot_wireframe(X,Y,F1)
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_zlabel(r"$f_2$")
    plt.xlabel(r"$\log_{10}(\tau_l=1./S_p)/\tau_s^{1}$")
    plt.ylabel(r"$\log_{10}(\tau_h= F_n/C_q)/\tau_s^{1}$")
    ax.scatter(X.ravel(),Y.ravel(),data_inf["f2"] )
    ax.plot_wireframe(X,Y,F2)
    plt.show()


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_zlabel(r"$f_1+f_2$")
    plt.xlabel(r"$\log_{10}(\tau_l=1./S_p)/\tau_s^{1}$")
    plt.ylabel(r"$\log_{10}(\tau_h= F_n/C_q)/\tau_s^{1}$")
    ax.scatter(X.ravel(),Y.ravel(),data_inf["f2"] +data_inf["f1"])
    ax.plot_wireframe(X,Y,F1+F2)
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_zlabel("Bound")
    plt.xlabel(r"$\log_{10}(\tau_l=1./S_p)/\tau_s^{1}$")
    plt.ylabel(r"$\log_{10}(\tau_h= F_n/C_q)/\tau_s^{1}$")
    ax.scatter(X.ravel(),Y.ravel(),data_inf["bound"] )
    ax.plot_wireframe(X,Y,B)
    plt.show()


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_zlabel("Finders")
    plt.xlabel(r"$\log_{10}(\tau_l=1./S_p)/\tau_s^{1}$")
    plt.ylabel(r"$\log_{10}(\tau_h= F_n/C_q)/\tau_s^{1}$")
    ax.scatter(X.ravel(),Y.ravel(),data_inf["find"] )
    ax.plot_wireframe(X,Y,B)
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_zlabel(r"$\tau_s^R$")
    plt.xlabel(r"$\log_{10}(\tau_l=1./S_p)/\tau_s^{1}$")
    plt.ylabel(r"$\log_{10}(\tau_h= F_n/C_q)/\tau_s^{1}$")
    ax.scatter(X.ravel(),Y.ravel(),data_inf["\\tau_s^R"] )
    ax.plot_wireframe(X,Y,T )
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

    #TS versus theory
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$S_n$")
    ax.set_zlabel(r"$\tau_s$")
    TSth=np.zeros((X.shape[0],X.shape[1]))
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            TSth[i,j]=tausr(*get_constants(PS_n=Y[i,j]) )
    ax.scatter(X.ravel(),Y.ravel(),TS[:,:].ravel())
    #ax.scatter(X.ravel(),Y.ravel(),F1th[:,:].ravel() ,color='r' )
    #ax.plot_wireframe(X,Y,F1th,color='r')
    ax.plot_wireframe(X,Y,TSth,color='g')
    plt.show()
    

    #f1 versus theory
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$S_n$")
    ax.set_zlabel(r"$f1$")
    F1th=np.zeros((X.shape[0],X.shape[1]))
    F1th2=np.zeros((X.shape[0],X.shape[1]))
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            F1th[i,j]=f1theo(X[i,j],*get_constants(PS_n=Y[i,j]) )
            F1th2[i,j]=fsbtheo(X[i,j],*get_constants(PS_n=Y[i,j]), taus=TS[i,j] )[0]
    ax.scatter(X.ravel(),Y.ravel(),F1[:,:].ravel())
    #ax.scatter(X.ravel(),Y.ravel(),F1th[:,:].ravel() ,color='r' )
    #ax.plot_wireframe(X,Y,F1th,color='r')
    ax.plot_wireframe(X,Y,F1th2,color='g')
    plt.show()
    
    #f2 versus theory
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$S_n$")
    ax.set_zlabel(r"$f2$")
    #ax.set_zlim3d(0, .2)
    F2th=np.zeros((X.shape[0],X.shape[1]))
    F2th2=np.zeros((X.shape[0],X.shape[1]))
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            F2th[i,j]=f2theo(X[i,j],*get_constants(PS_n=Y[i,j]) )
            F2th2[i,j]=fsbtheo(X[i,j],*get_constants(PS_n=Y[i,j]), taus=TS[i,j]  )[1]
    ax.scatter(X.ravel(),Y.ravel(),F2[:,:].ravel())
    #ax.scatter(X.ravel(),Y.ravel(),F2th[:,:].ravel() ,color='r' )
    #ax.plot_wireframe(X,Y,F2th,color='r')
    ax.plot_wireframe(X,Y,F2th2,color='g')
    plt.show()
    
    #bound versus theory
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$S_n$")
    ax.set_zlabel("Bound")
    #ax.set_zlim3d(0, .2)
    B=np.zeros((X.shape[0],X.shape[1]))
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            B[i,j]=fsbtheo(X[i,j],*get_constants(PS_n=Y[i,j]), taus=TS[i,j]  )[3]
    ax.scatter(X.ravel(),Y.ravel(),data["bound"].ravel())
    ax.plot_wireframe(X,Y,B,color='g')
    plt.show()
    
    #finders versus theory
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$S_n$")
    ax.set_zlabel("Finders")
    #ax.set_zlim3d(0, .2)
    FS=np.zeros((X.shape[0],X.shape[1]))
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            FS[i,j]=fsbtheo(X[i,j],*get_constants(PS_n=Y[i,j]), taus=TS[i,j]  )[2]
    ax.scatter(X.ravel(),Y.ravel(),data["find"].ravel())
    ax.plot_wireframe(X,Y,FS,color='g')
    plt.show()

    #Total fishers in school versus theory
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$S_n$")
    ax.set_zlabel("Bound fishers in school")
    #ax.set_zlim3d(0, .2)
    BHAT=np.zeros((X.shape[0],X.shape[1]))
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            BHAT[i,j]=fsbtheo(X[i,j],*get_constants(PS_n=Y[i,j]), taus=TS[i,j]  )[-1]
    ax.scatter(X.ravel(),Y.ravel(),data["otot"].ravel()-data["find"].ravel()*constants['PC_n'])
    ax.plot_wireframe(X,Y,BHAT,color='g')
    plt.show()


    #f1+f2 versus theory
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$S_n$")
    ax.set_zlabel(r"$f_1+f_2$")
    Fth=np.zeros((X.shape[0],X.shape[1]))
    Fth2=np.zeros((X.shape[0],X.shape[1]))
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            Fth[i,j]=f2theo(X[i,j],*get_constants(PS_n=Y[i,j]) )+f1theo(X[i,j],*get_constants(PS_n=Y[i,j]) )
            Fth2[i,j]=fsbtheo(X[i,j],*get_constants(PS_n=Y[i,j]), taus=TS[i,j]  )[-2]
    ax.scatter(X.ravel(),Y.ravel(),F1[:,:].ravel()+F2[:,:].ravel())
    #ax.scatter(X.ravel(),Y.ravel(),Fth[:,:].ravel() ,color='r' )
    ax.scatter(X.ravel(),Y.ravel(),Fth2[:,:].ravel() ,color='g' )
    #ax.plot_wireframe(X,Y,Fth,color='r')
    ax.plot_wireframe(X,Y,Fth2,color='g')
    plt.show()

    
    #Catch versus theory
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$S_n$")
    ax.set_zlabel(r"$H$")
    Hth=np.zeros((X.shape[0],X.shape[1]))
    Hth2=np.zeros((X.shape[0],X.shape[1]))
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            Hth[i,j]=Htheo(X[i,j],*get_constants(PS_n=Y[i,j]) )
            Hth2[i,j]=fsbtheo(X[i,j],*get_constants(PS_n=Y[i,j]), taus=TS[i,j]  )[-2]*constants["PC_q"]
            #Hth2[i,j]=(1-(1-B[i,j] - FS[i,j]) - B[i,j]* )*constants["PC_q"]
    ax.scatter(X.ravel(),Y.ravel(),H[:,:].ravel())
    #ax.scatter(X.ravel(),Y.ravel(),Hth[:,:].ravel(),color='r')
    #ax.plot_wireframe(X,Y,Hth,color='r')
    ax.scatter(X.ravel(),Y.ravel(),Hth2.ravel(),color='g')
    ax.plot_wireframe(X,Y,Hth2,color='g')
    plt.show()


    #Catch versus theory: Section
    fig = plt.figure()
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$H$")
    plt.scatter(X[:-1:2,-1],.5*(H[1::2,-1]+H[0:-1:2,-1]) )
    plt.plot(X[:,-1],Hth2[:,-1],color='g')
    plt.show()


if fig4opt_cn:
    #=========== LAMBDA AND CN ==================
    filename="Data_Fig4opt_cn"
    set_constants(filename)
    data=load_data(filename) 
    X=data['PC_lambda']
    Y=data['PC_n']
    TS=data['\\tau_s^R']
    F1=data['f1']
    F2=data['f2']
    H=data['Hrate']

    print constants["PS_n"]
    #f1 versus theory
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$C_n$")
    ax.set_zlabel(r"$f1$")
    F1th=np.zeros((X.shape[0],X.shape[1]))
    F1th2=np.zeros((X.shape[0],X.shape[1]))
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            F1th[i,j]=f1theo(X[i,j],*get_constants(PC_n=Y[i,j]) )
            F1th2[i,j]=fsbtheo(X[i,j],*get_constants(PC_n=Y[i,j]), taus=TS[i,j]  )[0]
    ax.scatter(X.ravel(),Y.ravel(),(F1.ravel()))
    #ax.scatter(X.ravel(),Y.ravel(),F1th[:,:].ravel() ,color='r' )
    #ax.plot_wireframe(X,Y,F1th,color='r')
    ax.plot_wireframe(X,Y,(F1th2),color='g')
    plt.show()
    
    #f2 versus theory
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$C_n$")
    ax.set_zlabel(r"$f2$")
    #ax.set_zlim3d(0, .2)
    F2th=np.zeros((X.shape[0],X.shape[1]))
    F2th2=np.zeros((X.shape[0],X.shape[1]))
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            F2th[i,j]=f2theo(X[i,j],*get_constants(PC_n=Y[i,j]) )
            F2th2[i,j]=fsbtheo(X[i,j],*get_constants(PC_n=Y[i,j]), taus=TS[i,j]  )[1]
    ax.scatter(X.ravel(),Y.ravel(),F2[:,:].ravel())
    #ax.scatter(X.ravel(),Y.ravel(),F2th[:,:].ravel() ,color='r' )
    #ax.plot_wireframe(X,Y,F2th,color='r')
    ax.plot_wireframe(X,Y,F2th2,color='g')
    plt.show()

    #bound versus theory
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$C_n$")
    ax.set_zlabel("Bound")
    #ax.set_zlim3d(0, .2)
    B=np.zeros((X.shape[0],X.shape[1]))
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            B[i,j]=fsbtheo(X[i,j],*get_constants(PC_n=Y[i,j]), taus=TS[i,j]  )[3]
    ax.scatter(X.ravel(),Y.ravel(),data["bound"].ravel())
    ax.plot_wireframe(X,Y,B,color='g')
    plt.show()
    
    #find versus theory
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$C_n$")
    ax.set_zlabel("Finders")
    #ax.set_zlim3d(0, .2)
    B=np.zeros((X.shape[0],X.shape[1]))
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            B[i,j]=fsbtheo(X[i,j],*get_constants(PC_n=Y[i,j]), taus=TS[i,j]  )[2]
    ax.scatter(X.ravel(),Y.ravel(),data["find"].ravel())
    ax.plot_wireframe(X,Y,B,color='g')
    plt.show()
    
    #f1+f2 versus theory
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$C_n$")
    ax.set_zlabel(r"$f_1+f_2$")
    Fth=np.zeros((X.shape[0],X.shape[1]))
    Fth2=np.zeros((X.shape[0],X.shape[1]))
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            Fth[i,j]=f2theo(X[i,j],*get_constants(PC_n=Y[i,j]) )+f1theo(X[i,j],*get_constants(PC_n=Y[i,j]) )
            Fth2[i,j]=fsbtheo(X[i,j],*get_constants(PC_n=Y[i,j]), taus=TS[i,j]  )[-2]
    ax.scatter(X.ravel(),Y.ravel(),F1[:,:].ravel()+F2[:,:].ravel())
    #ax.scatter(X.ravel(),Y.ravel(),Fth[:,:].ravel() ,color='r' )
    #ax.plot_wireframe(X,Y,Fth,color='r')
    ax.plot_wireframe(X,Y,Fth2,color='g')
    plt.show()
    

    
    #Catch versus theory
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$C_n$")
    ax.set_zlabel(r"$H$")
    Hth=np.zeros((X.shape[0],X.shape[1]))
    Hth2=np.zeros((X.shape[0],X.shape[1]))
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            Hth[i,j]=Htheo(X[i,j],*get_constants(PC_n=Y[i,j]) )
            Hth2[i,j]=fsbtheo(X[i,j],*get_constants(PC_n=Y[i,j]), taus=TS[i,j]  )[-2]*constants["PC_q"]
    ax.scatter(X.ravel(),Y.ravel(),H[:,:].ravel())
    #ax.scatter(X.ravel(),Y.ravel(),Hth[:,:].ravel(),color='r')
    #ax.plot_wireframe(X,Y,Hth,color='r')
    ax.scatter(X.ravel(),Y.ravel(),Hth2[:,:].ravel(),color='g')
    ax.plot_wireframe(X,Y,Hth2,color='g')
    plt.show()


