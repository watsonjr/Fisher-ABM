import numpy as np
import scipy.special as SS
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import *


#### SWITCHES FOR WHICH FIGURES TO PLOT ####

fig2a=1
fig2b=0
fig3=0


#=========== Extraction of constants from julia code ==================

fcst=open("../Constants.jl","r")
constants={}
for l in fcst:
	line=l.split()
	if line and ("const" in line[0] or "global" in line[0]):
		try:
			constants[line[1]]=eval(line[3].replace(';',''))
		except Exception as e:
			print e

GRD_mx = constants["GRD_nx"]*constants["GRD_dx"]


def get_constants(**kwargs):
	args=[ kwargs[i] if i in kwargs else constants[i] for i in 
		("PC_rp","PC_f","PC_q","PF_sig" ,"PF_n" ,"GRD_nx","GRD_dx","PC_v","PS_p")]
	return args

#=========== Analytical expressions ==================

def tausr_base(PC_rp,PC_f,PC_q,PF_sig ,PF_n,GRD_nx,GRD_dx,PC_v,PS_p):
	v=PC_v
	a=(PC_f+PF_sig) #min(PF_sig,PF_n*PC_f/(2*pi) )  #This correction applies for explicit fish
	b=GRD_nx*GRD_dx/2.
	t2=1./(1.-PC_rp) #avg time of straight flight (CHECKED NUMERICALLY)
	t1=1./(PC_rp) #avg time of wait (CHECKED NUMERICALLY)
	t2opt= (a/v * sqrt( log(b/a) -1./2.))
	popt=1/t2opt
	return v,a,b,t1,t2
	
def p_opt(PC_rp,PC_f,PC_q,PF_sig ,PF_n,GRD_nx,GRD_dx,PC_v,PS_p):
	v=PC_v
	a=PC_f+min(PF_sig,PF_n*PC_f/(2*pi) )
	b=GRD_nx*GRD_dx/2.
	t2opt= (a/v * sqrt( log(b/a) -1./2.))
	#if 1./t2opt>1:
	#	print (a,b,v,t2opt)
	popt=min(.98,1./t2opt)
	return popt

def tausr1(*args):  #Static mode
	PC_rp,PC_f,PC_q,PF_sig ,PF_n ,GRD_nx,GRD_dx,PC_v,PS_p=args
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
	
def tausr1b(*args):  #Static mode, infinite k
	PC_rp,PC_f,PC_q,PF_sig ,PF_n ,GRD_nx,GRD_dx,PC_v,PS_p=args
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
	PS_p=args[-1]
	t0= 1./PS_p #average time between jumps
	#return 1 + (exp(-t0/tt)-1)/(t0/tt)
	return (1- log( 1 + t0/tt) / (t0/tt))







#=========== FIGURES ==================

if fig2a:
    #=========== FIGURE 2A ==================

    ## Load in data
    TS = np.load("../Data/Data_firstpass.npy")
    xs = np.load("../Data/Data_firstpass_xs.npy")
    #print(TS)
    plt.plot(xs,TS)
    plt.plot(xs,[tausr(*get_constants(PC_rp=x)) for x in xs])
    plt.show()



if fig2a:
    #=========== FIGURE 2A ==================

    ## Load in data
    TS = np.load("../Data/Data_Fig2a.npy")
    F1 = np.load("../Data/Data_Fig2a_f1.npy")
    xs = np.load("../Data/Data_Fig2a_xs.npy")
    #print(TS)
    plt.plot(xs,TS)
    plt.plot(xs,[tausr(*get_constants(PC_rp=x)) for x in xs])
    plt.show()
    plt.plot(xs,F1)
    plt.plot(xs,[f1r(*get_constants(PC_rp=x)) for x in xs])
    plt.show()



if fig2b:
    #=========== FIGURE 2B ==================

    TS = np.load("../Data/Data_Fig2b.npy")
    xy = np.load("../Data/Data_Fig2b_xs.npy")

    #3dscatter
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel("F_sigma")
    plt.ylabel("C_f")
    #plt.zlabel("tau_s")
    ax.set_zlim(bottom=0, top=500)
    X,Y=xy[:,:,0].ravel(),xy[:,:,1].ravel()
    ax.scatter(X/GRD_mx,Y/GRD_mx,TS[:,:,0].ravel(),c=TS[:,:,0].ravel())

    #Wireframe
    X,Y=xy[:,:,0],xy[:,:,1]
    Z=X.copy()
    T=X.copy()
    for i in range(Z.shape[0]):
	    for j in range(Z.shape[1]):
		    Z[i,j]=TS[i,j,0]
		    args=get_constants(PF_sig=X[i,j],PC_f=Y[i,j])
		    args[0]=p_opt(*args) #optimal turn probability
		    T[i,j]=tausr(*args)

    ax.plot_wireframe(X/GRD_mx,Y/GRD_mx,T)
    plt.show()

    #Heatmap
    fig = plt.figure()
    plt.xlabel("F_sigma")
    plt.ylabel("C_f")
    ax = fig.add_subplot(111)
    ax.pcolor(X/GRD_mx, Y/GRD_mx, Z,  vmin=abs(Z).min(), vmax=500)#abs(Z).max())
    plt.show()

if fig3:
    #=========== FIGURE 3 ==================

    TS = np.load("../Data/Data_Fig3_noinf.npy")
    TSinf = np.load("../Data/Data_Fig3_inf.npy")
    xy = np.load("../Data/Data_Fig3_xs.npy")

    #3dscatter
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel("S_p")
    plt.ylabel("C_q")
    #plt.zlabel("tau_s")
    ax.set_zlim(bottom=0, top=1000)
    X,Y=xy[:,:,0].ravel(),xy[:,:,1].ravel()
    Z1,Z2=TS[:,:,0].ravel(),TSinf[:,:,0].ravel()
    X=1./X
    Y=constants["PF_n"]/Y
    X,Y=np.log10(X),np.log10(Y)
    ax.scatter(X,Y,Z1,c=Z1)
    ax.scatter(X,Y,Z2,c=Z2)


    #Wireframe
    X,Y=xy[:,:,0],xy[:,:,1]
    Z1=X.copy()
    Z2=X.copy()
    T=X.copy()
    for i in range(Z1.shape[0]):
	    for j in range(Z1.shape[1]):
		    Z1[i,j]=TS[i,j,0]
		    Z2[i,j]=TSinf[i,j,0]
		    args=get_constants(PS_p=X[i,j],PF_n=Y[i,j])
		    args[0]=p_opt(*args) #optimal turn probability
		    T[i,j]=tausr(*args)
    X=1./X
    Y=constants["PF_n"]/Y
    X,Y=np.log10(X),np.log10(Y)
    ax.plot_wireframe(X,Y,T)
    plt.show()

    #Heatmap
    fig = plt.figure()
    plt.xlabel("1./S_p")
    plt.ylabel("F_n/C_q")
    ax = fig.add_subplot(111)
    ax.pcolor(X, Y, Z1,  vmin=abs(Z1).min(), vmax=1000)#abs(Z).max())
    plt.show()
