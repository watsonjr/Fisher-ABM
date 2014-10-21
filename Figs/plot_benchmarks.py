import numpy as np
import scipy.special as SS
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import *
import os
from datatools import *


path="../Data/"

#### SWITCHES FOR WHICH FIGURES TO PLOT ####

firstpass=1
firstpass_ns=1
fig2a=1
fig2b=1
fig3=1
fig3bis=1 #Other quantities related to fig3
fig4opt=1
#=========== Extraction of constants from julia code ==================

if 0:
    ####OLD
    fcst=open("../Constants.jl","r")
    constants={}
    for l in fcst:
        line=l.split()
        if len(line)>1 and "=" in line[1]:
            try:
                constants[line[0]]=eval(line[2].replace(';',''))
            except Exception as e:
                print e
    GRD_mx = constants["GRD_nx"]*constants["GRD_dx"]
    fcst.close()
    
constants={}
def set_constants(filename):
    if not ".dat" in filename:
        filename=filename+".dat"
    fcst=open(os.path.join(path,filename ),"r")
    #constants={}
    global constants
    for l in fcst:
        line=l.split()
        if line and '#' in line[0]:
            continue
        if len(line)>1:
            exec("{} = {}".format(line[0],line[1]))
            exec("constants['{}'] = {}".format(line[0],line[1]))
    fcst.close()

def load_data(filename):
    if not ".npy" in filename:
        filename=filename+".npy"
    data= np.load(os.path.join(path,filename ))
    #get headers
    headers=[]
    fcst=open(os.path.join(path,filename.replace("npy","dat" )),"r")
    for l in fcst:
        line=l.split()
        if line and '#' in line[0]:
            line=[x.strip().replace('#','') for x in line ]
            headers=[ x for x in line if x]
            break
    datadict={}
    dim=len(data.shape)
    data=np.transpose(data,[dim-1]+range(dim)[:-1] ) #transpose to get dataset as first index
    for i,h in enumerate(headers):
        datadict[h]=data[i]
    return datadict

set_constants("Data_params")

def get_constants(**kwargs):
    args=[ kwargs[i] if i in kwargs else constants[i] for i in 
        ("PC_rp","PC_f","PC_q","PF_sig" ,"PF_n" ,"GRD_nx","GRD_dx","PC_v","PS_p","PS_n")]
    return args

#=========== Analytical expressions ==================

def domain_size(PC_rp,PC_f,PC_q,PF_sig ,PF_n,GRD_nx,GRD_dx,PC_v,PS_p,PS_n):
    b=((PF_sig+PC_f)/2+GRD_nx*GRD_dx/2./sqrt(PS_n) )
    return b

def tausr_base(PC_rp,PC_f,PC_q,PF_sig ,PF_n,GRD_nx,GRD_dx,PC_v,PS_p,PS_n):
    v=PC_v
    a=(PC_f+PF_sig) #min(PF_sig,PF_n*PC_f/(2*pi) )  #This correction applies for explicit fish
    b=domain_size(PC_rp,PC_f,PC_q,PF_sig ,PF_n,GRD_nx,GRD_dx,PC_v,PS_p,PS_n)
    PC_rp=max(PC_rp,0.01)
    t2=1./(1.-PC_rp) #avg time of straight flight (CHECKED NUMERICALLY)
    t1=1./(PC_rp) #avg time of wait (CHECKED NUMERICALLY)
    if log(b/a)>.5:
        t2opt= (a/v * sqrt( log(b/a)-1./2.))
    else:
        t2opt=0.0001
    popt=1/t2opt
    return v,a,b,t1,t2
    
def p_opt(PC_rp,PC_f,PC_q,PF_sig ,PF_n,GRD_nx,GRD_dx,PC_v,PS_p,PS_n):
    v=PC_v
    a=PC_f+PF_sig#min(PF_sig,PF_n*PC_f/(2*pi) )
    b=domain_size(PC_rp,PC_f,PC_q,PF_sig ,PF_n,GRD_nx,GRD_dx,PC_v,PS_p,PS_n)
    if log(b/a)>.5:
        t2opt= (a/v * sqrt( log(b/a)-1./2.))
    else:
        t2opt=0.0001
    #if 1./t2opt>1:
    #    print (a,b,v,t2opt)
    popt=min(.98,1./t2opt)
    return popt

def tausr1(*args):  #Static mode
    PC_rp,PC_f,PC_q,PF_sig ,PF_n ,GRD_nx,GRD_dx,PC_v,PS_p,PS_n=args
    v,a,b,t1,t2=tausr_base(*args)
    k=1./(1.-PC_q) #rate of catching, to specify better
    xy=sqrt( 2*k*t1/(1+k*t1) )/v/t2
    x=a*xy
    y=b*xy
    result1=xy**1/a*(1+k*t1)*(b**2-a**2)**2*SS.iv(0,x)/SS.iv(1,x)
    result2=  8*b**2 + xy**2 * (1+k*t1)*(4*b**4*log(b/a) + (b**2-a**2)*(a**2-3*b**2+8/xy**2))
    result3= result1+ result2/4.
    result4= (t1+t2)/(2*k*t1*b**2) * result3
    return result4
    
def tausr1b(*args):  #Static mode, infinite k (USE THIS ONE!!!!~)
    PC_rp,PC_f,PC_q,PF_sig ,PF_n ,GRD_nx,GRD_dx,PC_v,PS_p,PS_n=args
    v,a,b,t1,t2=tausr_base(*args)
    xy=sqrt(2)/v/t2
    x=a*xy
    y=b*xy
    result1=xy**1/a*t1*(b**2-a**2)**2*SS.iv(0,x)/SS.iv(1,x)
    result2=xy**2 *t1*(4*b**4*log(b/a) + (b**2-a**2)*(a**2-3*b**2+8/xy**2))
    result3= result1 + result2/4.
    result4= (t1+t2)/(2*t1*b**2) * result3
    return result4
    
    
def tausr2(*args):  #Diffusive mode
    v,a,b,t1,t2=tausr_base(*args)
    D=1.
    Db=v**2*t2/2
    alp=sqrt(1./(D*t1) + 1./(Db*t2) )
    elem1=SS.iv(1,b*alp) * SS.kv(1,a*alp)-SS.iv(1,a*alp) * SS.kv(1,b*alp)
    elem2=SS.iv(1,b*alp) * SS.kv(0,a*alp)+SS.iv(0,a*alp) * SS.kv(1,b*alp)
    Lelem1=SS.iv(0,a/sqrt(Db*t2) )* elem1
    Lelem2= alp * sqrt(Db*t2) *SS.iv(1,a/sqrt(Db*t2) )* elem2
    Lp=Lelem1 + Lelem2
    Lm= Lelem1-Lelem2
    M= SS.iv(0,a/sqrt(Db*t2) )* elem2  -4 * a**2* sqrt(Db*t2) /(alp* (b**2-a**2)**2) *SS.iv(1,a/sqrt(Db*t2) )* elem1
    result1=a*alp*(b**2/a**2-1) * M/2/Lp - Lm / Lp
    result2= (3-4*log(b/a))*b**4 - 4*a**2 *b**2 + a**4
    result3= result1 - a**2 * D *t1/ (8*Db*t2) * result2 / (b**2-a**2)
    result4=(t1+t2) * (1-a**2/b**2)/(alp**2 * D *t1)**2 * result3
    return result4

def tausr3(*args):  #Minimal approximation
    v,a,b,t1,t2=tausr_base(*args)
    x= a * sqrt(2.)/v/t2
    result1= (b**2-a**2)**2/ (sqrt(2) *v*a*b**2 ) *  SS.iv(0,x)/SS.iv(1,x)
    result2= b**4 * log(b/a) + (b**2-a**2) * (a**2 - 3 *b**2 + 4 * v**2 *t2**2)/4.
    return   result1 + result2/(v**2 * t2 * b**2)

def tausr(*args):
    tt= tausr1b(*args)
    return tt
    

def f1r(*args):
    tt= tausr1b(*args)
    PC_rp,PC_f,PC_q,PF_sig ,PF_n ,GRD_nx,GRD_dx,PC_v,PS_p,PS_n=args
    th= PF_n/PC_q
    if PS_p>0:
        tl= 1./PS_p #average time between jumps
        if th<tl:
            print("quick harvest")
    else:
        tl= 10000000000
    t0=min(tl,th)
    #return 1 + (exp(-t0/tt)-1)/(t0/tt)
    return (1- log( 1 + t0/tt) / (t0/tt))







#=========== FIGURES ==================

if firstpass:
    filename="Data_firstpass"
    set_constants(filename)
    data=load_data(filename) 
    TS=data['\\tau_s^R']
    xs=data['PC_rp']
    #print(TS)
    plot(TS,[tausr(*get_constants(PC_rp=x)) for x in xs],xs=xs)


if firstpass_ns:
    filename="Data_firstpass_ns"
    set_constants(filename)
    data=load_data(filename) 
    TS=data['\\tau_s^R']
    dist=data['dist']
    xs=data['PS_n']
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



if fig2b:
    #=========== FIGURE 2B ==================

    filename="Data_Fig2b"
    set_constants(filename)
    data=load_data(filename) 
    TS=data['\\tau_s^R']
    F1=data["f1"]
    X=data['PF_sig']
    Y=data['PC_f']

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
            args=get_constants(PF_sig=X[i,j],PC_f=Y[i,j])
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
   # X,Y=np.log10(X),np.log10(Y)

    taus1=tausr(*get_constants(**constants))
#    ax.set_zlim(bottom=-1, top=1)
    ax.scatter(X.ravel()/taus1,Y.ravel()/taus1,Z1/Z2-1,c=Z1/Z2-1)
    ax.set_zlabel(r"$VOI= \tau_s^R/\tau_s^I-1$")
    plt.show()
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\log_{10}(\tau_l=1./S_p)/\tau_s^{1}$")
    plt.ylabel(r"$\log_{10}(\tau_h= F_n/C_q)/\tau_s^{1}$")
    ax.set_zlim(bottom=0, top=1000)
    ax.scatter(X.ravel()/taus1,Y.ravel()/taus1,Z1,c='r')
    ax.scatter(X.ravel()/taus1,Y.ravel()/taus1,Z2,c='b')
    ax.set_zlabel(r"$\tau_s$")

    #Wireframe
    T=X.copy()
    for i in range(T.shape[0]):
        for j in range(T.shape[1]):
            args=get_constants(PS_p=data_inf['PS_p'][i,j],PC_q=data_inf['PC_q'][i,j] )
            args[0]=p_opt(*args) #optimal turn probability
            T[i,j]=tausr(*args)

    ax.plot_wireframe(X/taus1,Y/taus1,T)
    plt.show()

    #Heatmap
    #Z1=TS
    #Z2=TSinf
    #fig = plt.figure()
    #plt.xlabel(r"$\log_{10}\tau_l=1./S_p$")
    #plt.ylabel(r"$\log_{10}\tau_h= F_n/C_q$")
    #ax = fig.add_subplot(111)
    #ax.pcolor(X, Y, Z1-Z2,  vmin=(Z1-Z2).min(), vmax=(Z1-Z2).max())#abs(Z).max())
    #plt.show()
    
    
if fig3bis:
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
    Catchinf=data_inf["H"]
    X=data_inf['PS_p']
    Y=data_inf['PC_q']

    X=1./X
    Y=constants["PF_n"]/Y
    X,Y=np.log10(X),np.log10(Y)
    
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
    
    
    #CATCH
    taus1=tausr(*get_constants())

    T1,T2=Catch[:,:].ravel(),Catchinf[:,:].ravel()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel(r"$\log_{10}(\tau_l=1./S_p)/\tau_s^{1}$")
    plt.ylabel(r"$\log_{10}(\tau_h= F_n/C_q)/\tau_s^{1}$")
    Z=T2/T1-1
    ax.scatter(X.ravel()/taus1,Y.ravel()/taus1,Z,c=Z)
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
    H=data['Hdist']

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
    
    mycmap = plt.cm.get_cmap('Greys')
    #mycmap.set_under('w')
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim([1,10])
    plt.ylim([1,20])
    plt.pcolor(X,Y,TS[:,:], cmap=mycmap,vmin =min(TS[:,:].ravel()), vmax=max(TS[:,:].ravel()))
    plt.colorbar()
    plt.show()
    
    
    filename="Data_Fig4opt_cliq"
    set_constants(filename)
    data=load_data(filename) 
    X=data['PC_ncliq']
    Y=data['PS_n']
    TS=data['\\tau_s^R']
    H=data['Hdist']

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

    mycmap = plt.cm.get_cmap('Greys')
    #mycmap.set_under('w')
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim([1,10])
    plt.ylim([1,20])
    plt.pcolor(X,Y,TS[:,:], cmap=mycmap,vmin =min(TS[:,:].ravel()), vmax=max(TS[:,:].ravel()))
    plt.colorbar()
    plt.show()
