from plot_core import *

#=========== Analytical expressions ==================

def domain_size(PC_rp,PC_f,PC_q,PF_sig ,PF_n,GRD_nx,GRD_dx,PC_v,PS_p,PS_n,PC_n):
    X=GRD_nx*GRD_dx
    b=(PF_sig+PC_f+X/sqrt(PS_n ) )/2
    return b

def tausr_base(PC_rp,PC_f,PC_q,PF_sig ,PF_n,GRD_nx,GRD_dx,PC_v,PS_p,PS_n,PC_n):
    v=PC_v
    a=(PC_f+PF_sig) #min(PF_sig,PF_n*PC_f/(2*pi) )  #This correction applies for explicit fish
    b=domain_size(PC_rp,PC_f,PC_q,PF_sig ,PF_n,GRD_nx,GRD_dx,PC_v,PS_p,PS_n,PC_n)
    PC_rp=max(PC_rp,0.01)
    t2=1./(PC_rp) #avg time of straight flight (CHECKED NUMERICALLY)
    t1=1#/(1.-PC_rp) #avg time of wait (CHECKED NUMERICALLY)
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
    k=1./t1 /(1.-PC_q) #rate of catching, to specify better
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

def fsbtheo(PC_lambda,*args,**kwargs):
    args=list(args)
    PC_rp,PC_f,PC_q,PF_sig ,PF_n ,GRD_nx,GRD_dx,PC_v,PS_p,PS_n,PC_n=args
#    args[-2]=sqrt(PS_n)
    tl,th,t0=taucalc(*args)
    ts=kwargs.get('taus',tausr(*args))
    lam=PC_lambda**2
    #lam=1
    N=PC_n
    xi=t0/ts
    tl*=2
    if lam>0:
        #lam=0.5
        a,b,c,d = lam,1-lam,-N*(1-xi),-N**2*xi
        def optrad(f):
            #return a*f**3+b*f**2+c*f+d
            #s=f-N/(1.+th/ts)
            
            s=f/ts -.5*(N-f)/th
            s/=(1/ts + .5/tl + .5/th)
            if s<0:
                s=0
            return (N-f) -lam*s*(f-s)
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
    s=f/ts - (N-f)/th
    s/=(1/ts + .5/tl + .5/th)
    b=PC_n-f
    f,s,b=f/N,s/N,b/N
    #Def: f1+f2 = 1-f0
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
    fnc=tnc/(t0+tnc)
    f0=(f-s)+b*fnc
    
    def get_f(s):
        return (N + lam * s**2)/(1+lam*s)
    def get_b(s):
        return N-get_f(s)
    def get_no(s):
        return PS_n * (1-exp(-s/PS_n) )
    def get_bhat(s):
        f=get_f(s)
        #term=s/N*sqrt(PC_rp/2*get_no(s))
        #term=2*PC_v/(GRD_nx*GRD_dx)
        #bhat=(N-f)/(1+(f-s)/(ts*s)/term   )
        bhat= th*get_no(s)*( (N-s)/ (ts*s * (1+lam*s) ) - 1./tl) -s
        
        chi=sqrt(PC_rp/2/get_no(s))*s/N
        t2=N-f-s-(f-s)/ts/chi
        t3=s*(N-f)
        #bhat=-t2/2 + sqrt(t2**2/4 - t3)
        #bhat=chi*s*(N-f)/(chi*s + (f-s)/ts/(get_b(s)/s) )
        return bhat
        
    if 0:
        s*=N

        f=get_f(s)
        b=get_b(s)
        no=get_no(s)
        bhat=get_bhat(s)
        omean=(s+bhat)/no
        f0=1.-(s+bhat)*1./N
        #print omean,no,s
        f1=s/N*exp(-bhat/s)*exp(-s/PS_n) # no/PS_n* omean* exp(-omean)/(1-exp(-omean))
        f2=1-f0-f1
        s/=N
        b/=N

    f,s,b=f*N,s*N,b*N
    return f1, f2,s,b,1-f0,get_bhat(s)
    
    
def fsbtheo_test(PC_lambda,*args):
    args=list(args)
    PC_rp,PC_f,PC_q,PF_sig ,PF_n ,GRD_nx,GRD_dx,PC_v,PS_p,PS_n,PC_n=args
    tl,th,t0=taucalc(*args)
    ts=tausr(*args)
    lam=PC_lambda**2
    N=PC_n
    D=2*PC_v**2/PC_rp
    def get_ns(vec):
        f,s,b,i=vec
        ns=(f/ts-s/th)/(f/(ts*PS_n)+1/tl)+s
        return ns
    def optrad(vec):
        #print xxx
        f,s,b,i=vec
        ns=get_ns(vec)
        d=(f+b)/N * GRD_nx*GRD_dx/2 + (i+s)/N*(PC_f+PF_sig) #typical distance between fishers
        fs=f*1/ts
        bs=b*PC_v/d
        si=s*(1./tl + s/ns / th)
        iff=i*D/PC_f**2
        fb=lam *f
        res= np.abs(np.array([fs+fb-iff,fs+bs-si,fb-bs,si-iff , f+b+s+i-N, f-abs(f),s-abs(s),b-abs(b),i-abs(i),ns-abs(ns)] ))
        #print res
        return res
    res = opt.leastsq(optrad,np.array( [N/2.,N/2.,N/4.,N/4. ]))#,method='TNC',bounds=[(0,N),(0,N),(0,N),(0,N)] )#, f_tol=1e-14)
    f,s,b,i=res[0]/N
    ns=get_ns(res[0])
    #print f,s,b,i,ns
    #Def: f1+f2 = 1-f0
    f0=f
    mu=s*N/ns
    if mu<0:
        mu=0
    f1=s*exp(-mu)
    f2=s-f1
    
    
    return f1, f2,s,b,1-f0
    
def fsbtheo(PC_lambda,*args,**kwargs):

    args=list(args)
    PC_rp,PC_f,PC_q,PF_sig ,PF_n ,GRD_nx,GRD_dx,PC_v,PS_p,PS_n,PC_n=args
    tl,th,t0=taucalc(*args)
    ts=kwargs.get('taus',tausr(*args))
    lam=PC_lambda**2
    N=PC_n
    xi=t0/ts
    tl*=2
    def get_f(s):
        return (N + lam * s**2)/(1+lam*s)
    def get_b(s):
        return N-get_f(s)
    def get_no(s):
        return PS_n * (1-exp(-s/PS_n) )
    def get_bhat(s):
        f=get_f(s)
        #term=s/N*sqrt(PC_rp/2*get_no(s))
        #term=2*PC_v/(GRD_nx*GRD_dx)
        #bhat=(N-f)/(1+(f-s)/(ts*s)/term   )
        bhat= th*get_no(s)*( (N-s)/ (ts*s * (1+lam*s) ) - 1./tl) -s
        return bhat
    def get_s(f):
        return f/2. + sqrt(f**2-4*(N-f)/lam)/2.

    def optrad(s):
        #print s
        if s<0:
            s=exp(s)/10
        bhat=get_bhat(s)
        f=get_f(s)
        no=get_no(s)

        #print '    ',bhat
        term1= (N-s)*bhat/ts
        term2= (N-bhat)*(1+lam*s)
        term2+= N+lam*s**2
        term2*=s
        term2*=2*PC_v/GRD_nx/GRD_dx*(s+bhat)/N#max((s+bhat)/N*sqrt(PC_rp/2*get_no(s) ),2*PC_v/GRD_nx/GRD_dx)
        
        term1=(f-s)*bhat/ts
        term2=2*PC_v/GRD_nx/GRD_dx * (N-f-bhat)*s
        
        #ALTERNATIVE
        #term1=(N-s)/ts
        b=(N-f)
        #b=bhat
        #print s,bhat,N
        #term2=s*(1+lam*s)*(1./tl+(s+bhat )/(th*no) )
        
        #term1=(f-s)/ts
        #term2=s*(1./tl+(s+bhat)/(th*no))
        
        return term2-term1
    s = opt.newton_krylov(optrad,N/ 4., f_tol=1e-14)
    #print s
    f=get_f(s)
    b=get_b(s)
    no=get_no(s)
    bhat=get_bhat(s)
    omean=(s+bhat)/no
    f0=1.-(s+bhat)*1./N
    #print omean,no,s
    f1=s/N*exp(-bhat/s)*exp(-s/PS_n)  #no/PS_n* omean* exp(-omean)/(1-exp(-omean))
    f2=1-f0-f1
    
    return f1, f2,s,b,1-f0,bhat*N
