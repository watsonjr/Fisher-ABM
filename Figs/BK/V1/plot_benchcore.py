import numpy as np
import scipy.special as SS
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import *
import os
from datatools import *
import cmath
import scipy.optimize as opt

path="../Data/"

#=========== Extraction of constants from julia code ==================


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
        try:
            datadict[h]=data[i]
        except:
            print "ERROR in load_data: ", headers, data.shape[0]
    return datadict

def get_constants(**kwargs):
    args=[ kwargs[i] if i in kwargs else constants[i] for i in 
        ("PC_rp","PC_f","PC_q","PF_sig" ,"PF_n" ,"GRD_nx","GRD_dx","PC_v","PS_p","PS_n","PC_n")]
    return args

#=========== Analytical expressions ==================

def domain_size(PC_rp,PC_f,PC_q,PF_sig ,PF_n,GRD_nx,GRD_dx,PC_v,PS_p,PS_n,PC_n):
    b=(PF_sig+PC_f+GRD_nx*GRD_dx/sqrt(PS_n) )/2
    return b

def tausr_base(PC_rp,PC_f,PC_q,PF_sig ,PF_n,GRD_nx,GRD_dx,PC_v,PS_p,PS_n,PC_n):
    v=PC_v
    a=(PC_f+PF_sig) #min(PF_sig,PF_n*PC_f/(2*pi) )  #This correction applies for explicit fish
    b=domain_size(PC_rp,PC_f,PC_q,PF_sig ,PF_n,GRD_nx,GRD_dx,PC_v,PS_p,PS_n,PC_n)
    PC_rp=max(PC_rp,0.01)
    t2=1./(PC_rp) #avg time of straight flight (CHECKED NUMERICALLY)
    t1=1./(1-PC_rp) #avg time of wait (CHECKED NUMERICALLY)
    return v,a,b,t1,t2
    
def p_opt(PC_rp,PC_f,PC_q,PF_sig ,PF_n,GRD_nx,GRD_dx,PC_v,PS_p,PS_n,PC_n):
    v=PC_v
    a=PC_f+PF_sig#min(PF_sig,PF_n*PC_f/(2*pi) )
    b=domain_size(PC_rp,PC_f,PC_q,PF_sig ,PF_n,GRD_nx,GRD_dx,PC_v,PS_p,PS_n,PC_n)
    if log(b/a)>.5:
        t2opt= (a/v * sqrt( log(b/a)-1./2.))
    else:
        t2opt=0.0001
    #if 1./t2opt>1:
    #    print (a,b,v,t2opt)
    popt=min(.98,1./t2opt)
    return popt

def tausr1(*args):  #Static mode
    PC_rp,PC_f,PC_q,PF_sig ,PF_n ,GRD_nx,GRD_dx,PC_v,PS_p,PS_n,PC_n=args
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
    PC_rp,PC_f,PC_q,PF_sig ,PF_n ,GRD_nx,GRD_dx,PC_v,PS_p,PS_n,PC_n=args
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
    Db=v**2*t1/2
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
    

def f1r(*args,**kwargs):
    tt= tausr1b(*args)
    PC_rp,PC_f,PC_q,PF_sig ,PF_n ,GRD_nx,GRD_dx,PC_v,PS_p,PS_n,PC_n=args
    tl,th,t0=taucalc(*args)
    th= kwargs.get('tauh',th)
    x=t0/tt
    return t0/(t0+tt)
    #return 1 + (exp(-x)-1)/(x)
    #return (1- log( 1 +x) / (x))


def f2rec(level,PC_lambda,*args):
    #RECURSIVE
    PC_rp,PC_f,PC_q,PF_sig ,PF_n ,GRD_nx,GRD_dx,PC_v,PS_p,PS_n,PC_n=args
    GRD_mx=GRD_nx*GRD_dx
    f1=f1r(*args)
    nu=0.
    if PC_n==1:
        lam=0
    else:
        lam=PC_lambda**(2 / (PC_n-1))
    chi=min(1,nu*PC_v/PS_p/GRD_mx)
    a=(PC_f+PF_sig)
    f2random= (1-lam)*f1**2 /PS_n 
    f2random= 1-(1-f2random)**(PC_n-1) #N Fishers extension
        
    if level<=0:
        #First approximation
        f2= f2random + f1*((1-lam)*f1 /PS_n +(1-nu)* lam  )
    else:
        f2=f2rec(level-1,PC_lambda,*args)
        #Time spent not collapsed
        tl,th,t0=taucalc(*args)
        if f2>0:
            tnc=t0*(1-f2)/f2
            D=2*PC_v**2/PC_rp
            fnc=1/(t0+tnc) * sqrt(D*tnc)/PC_v
        else:
            fnc=1
        
        zargs=list(args[:])
        zargs[1]=PC_f*(1+lam*(1-f2  ))
        f1bis=f1r(*zargs,tauh=PF_n/PC_q/(1+f2*(PC_n-1) ))
        
    #    f1bis=f1
        #f1=f1bis
        #f1=f1+lam*(1-lam)*2*f1
        f2info= (1-nu)* lam *f1bis*(1-fnc) #no difference wit number of fishers
        f2= f2random/2+f2info
    return f2

def f2theo(*args):
    return f2rec(2,*args)

def f1theo(PC_lambda,*args):
    PC_rp,PC_f,PC_q,PF_sig ,PF_n ,GRD_nx,GRD_dx,PC_v,PS_p,PS_n,PC_n=args
    f2=f2theo(PC_lambda,*args)
    if PC_n==1:
        lam=0
    else:
        lam=PC_lambda**(2 / (PC_n-1))
    tl,th,t0=taucalc(*args)
    f1=f1r(*args) 
    v,a,b,t1,t2=tausr_base(*args)
    GRD_mx=GRD_nx*GRD_dx
    if f2>0:
        tnc=t0*(1-f2)/f2
        D=2*PC_v**2/PC_rp
        fnc=1/(t0+tnc) * sqrt(D*tnc)/PC_v
    else:
        fnc=1
    fnc=(fnc)**(PC_n-1)
    return f1*((1-lam)*(1-f1 *(PC_n-1) /PS_n) + lam*fnc)# *max(0,1-v*GRD_mx/sqrt(PS_n)/2/t0))#*sqrt(PS_n)) )


def Htheo(PC_lambda,*args):
    PC_rp,PC_f,PC_q,PF_sig ,PF_n ,GRD_nx,GRD_dx,PC_v,PS_p,PS_n,PC_n=args
    f2=f2theo(PC_lambda,*args)
    f1=f1theo(PC_lambda,*args)#+f2/2
    C1=PC_q
    #tauh=PF_n/PC_q/PC_n    
    C2=PC_q#max(min(PC_q,PC_q*tauh*PS_p),PC_q/2) 
    return (C1*f1+ C2*f2)
    
    
def taucalc(*args):
    PC_rp,PC_f,PC_q,PF_sig ,PF_n ,GRD_nx,GRD_dx,PC_v,PS_p,PS_n,PC_n=args
    th= PF_n/PC_q
    if PS_p>0:
        tl= 1./PS_p #average time between jumps
    else:
        tl= 10000000000
    t0=min(tl,th)
    return tl, th, t0

def fsbtheo(PC_lambda,*args):
    args=list(args)
    PC_rp,PC_f,PC_q,PF_sig ,PF_n ,GRD_nx,GRD_dx,PC_v,PS_p,PS_n,PC_n=args
#    args[-2]=sqrt(PS_n)
    tl,th,t0=taucalc(*args)
    ts=tausr(*args)
    lam=PC_lambda**2
    #lam=1
    xi=ts/t0
    N=PC_n
    if lam>0:
        #lam=0.5
        a,b,c,d = lam,1-lam,-N*(1-xi),-N**2*xi
        def optrad(f):
            #return a*f**3+b*f**2+c*f+d
            #s=f**2/(N*xi+f)
            s=f-N/(1.+th/ts)
            #if th*s/(s+N*f)>tl:
             #   s=f/(1+tl/ts)
            return (N-f) -lam*s*(f-1)
        #a,b,c,d = 1+lam*(N*xi-1),N*(xi*(2-lam)-1),N**2*xi*(xi-2),-N**3*xi**2
        D0=b**2 - 3*a*c
        D1=2*b**3 - 9 *a*b*c + 27*a**2 * d
        D=cmath.sqrt(D1**2-4*D0**3)
        C=(.5*D1+.5*D ) **(1./3.)
        uk=1
        #f=-1/(3*a) *( b +uk* C + D0/C/uk)
        #f=[x for x in sorted(np.roots([a,b,c,d]),key=lambda e:abs(e.imag)) if float(x).real>0][0]
        f = opt.newton_krylov(optrad, N/2., f_tol=1e-14)
    else:
        f=N
    s=f**2/(N*xi+f)
    b=PC_n-f
    f,s,b=f/N,s/N,b/N
    #Def: f1+f2 = 1-f0
    f0=(f-s)
    #f1=s*exp(-s*N/PS_n)*exp(-b/s)
    if N>1:
        f1=s*(1-s/PS_n )**(N-1) *exp(-b/s)# ( max(0,1-1./s/N) )**(b*N)
    else:
        f1=s
    f2_rnd=(s-f1)
    #f2_inf=s*(1-exp(- b/f ))
    
    if N>1:
        f2_inf=b *s/f #*(1 - (1-b/( f*(N-1) ) )**(N-1) )
    else:
        f2_inf=0
    f2=f2_rnd + f2_inf
    
    if f2>0:
        tnc=t0*(1-f2)/f2
    else:
        tnc=t0*1000
    fnc=1/(t0+tnc) * tnc
    
    return f1, f2,s,b,1-(f-s)-b*fnc #f1+f2
