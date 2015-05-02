from datatools import *
from plot_benchcore import *
from random_partition import random_partition

ni,nj=range(1,9),range(1,9)

USE_RANDOM=0
DISCRETIZE_SIMU=0 #Discrete jumps

RANDOM_PARTITION=1 #Parition fishers in random cliques
nens= 50 #number of ensembles

filename="Data_anaflux"
set_constants(filename)
data=load_data(filename) 

keynames=("f","s","b","no","o * no","tauh" ,"taul")

evtdic={e:{} for  e in ("found_school","left_school","lost_school", "bound")}

constants={}
states={} #States throughout the simulation
series={} #Time series of te states
distdic={} #Typical distance between fishers

cliqH={} #Typical catch by clique size


def runavg(array,window):
    return np.array([np.mean(array[z:z+window]) for z in range(array.shape[0]/window)] )

def cumavg(array):
    #Cumsum rescaled by linear
    res=[]
    cur=0
    for i in range(array.shape[0]):
        cur*=i/(i+1.)
        cur+=array[i]/(i+1.)
        res.append(cur)
    return np.array(res)
    
for i in ni:
    print i
    for j in nj:
        d=np.load("../Data/dist-{}-{}.npy".format(i,j))#[::5]
        s=np.load("../Data/states-{}-{}.npy".format(i,j))#[::5] 
        cliq=np.load("../Data/cliq-{}-{}.npy".format(i,j))#[::5] 
        cliq[cliq==0]=np.NaN
        cliqH[(i,j)]=cliq
        
        series[(i,j)]=s
        
        for x in range(s.shape[0]):
            key=tuple(s[x] )
            key+=(d[x,1],i,j)
            states.setdefault(key,0)
            states[key]+=1.
            
        lines=open("../Data/turn{}-{}.dat".format(i,j),'r').readlines()
        nturns,ts,tstheo=[float(l.strip()) for l in lines[:3] ]
        constants[(i,j)]={l.split()[0]:eval(l.split()[1])  for l in lines[3:] }
        constants[(i,j)]["taus"]=ts
        constants[(i,j)]["taustheo"]=tstheo 
        
        #scatter(d[:,0],d[:,1]/constants[(i,j)]["GRD_mx2"]*sqrt(2.))
        distdic[(i,j)]=mhist(d[:,1]/constants[(i,j)]["GRD_mx2"]*sqrt(2.),normed=1,bins=100)
        xs,ys=distdic[(i,j)]
        #plot(xs,[exp(-(x-.8)**2/2*40) for x in xs],hold=1)
        #plot(xs,ys,log='y')
        #scatter(d[:,0],d[:,1]/constants[(i,j)]["GRD_mx2"]*sqrt(2.))
        
        #USELESS TEST (i wanted to relate typical distance d to tausr)
        def maketausr(sig):
            const={}
            const.update(constants[(i,j)])
            del const['PF_sig']
            del const['PS_n']
            return tausr(*get_constants( PF_sig=sig,PS_n=1,**const ))
        #res= [ 1./maketausr( x*constants[(i,j)]["GRD_mx2"]*sqrt(2.)) for x in xs]
        #res=array(res)
        #plot(xs[:-1],res[1:]-res[:-1],log='y')
                           
        for  e in ("found_school","left_school","lost_school", "bound"):
            localdic={}
            arr=np.load("../Data/{}-{}-{}.npy".format(e,i,j))

            for x in range(arr.shape[0]):
                key=tuple(arr[x])
                key+=(0,i,j )
                evtdic[e].setdefault(key,0.)
                evtdic[e][key]+=1./nturns

#Find which variables change along the axes
x_ax={}         
y_ax={}
for key,val0 in constants.values()[0].iteritems():
    key_array=np.zeros((len(ni),len(nj))+np.array(val0).shape) 
    for idxi,i in enumerate(ni):
        for idxj,j in enumerate(nj):
            key_array[idxi,idxj]=constants[(i,j)][key]
    if (key_array[:,0]==key_array[:,1]).all() and not (key_array[0,:]==key_array[1,:]).all():
        y_ax[key]=key_array
    if (key_array[0,:]==key_array[1,:]).all() and not (key_array[:,0]==key_array[:,1]).all():
        x_ax[key]=key_array
    
    
def get_stateproba(k,stateproba):
    #Get probability of a given fisher state among full time series
    pro=1.
    for col in range(2):
        x=k[col]
        xproba=bisct.bisect_left( stateproba[col][0],x )
        p_x=stateproba[col][1][xproba]/np.mean(stateproba[col][1])
        pro*=p_x
    return pro
    

def binner(dic, col=0,nbins=100):
    keys=np.array(dic.keys())
    krange=keys[:,col]
    if len(set(krange))<nbins:
        bins=sorted(set(krange))
    else:
        kmin,kmax=np.min(krange),np.max(krange)
        bins=np.linspace(kmin,kmax,num=nbins)
    binned={}
    for k in dic:
        xbin=bisct.bisect_left(bins,k[col]) 
        binned.setdefault(xbin,[]).append(k)
    return bins,binned
    
def binsumbase(dic,col=0,nbins=100):
    bins,binned=binner(dic,col,nbins)
    binval=np.zeros(len(bins))
    for xbin in binned:
        binval[xbin]=sum([dic[k] for k in binned[xbin]])
    return bins,binval


def binsum(dic, col=0,nbins=100,compare_to=None,stateproba=None):
    bins,binned=binner(dic,col,nbins)
    binval=np.zeros(len(bins))
    occur=np.zeros(len(bins))
    theor=np.zeros(len(bins))
    comp={}
    for xbin in binned:
        binval[xbin]=sum([dic[k]/get_stateproba(k,stateproba) for k in binned[xbin]])
        for k in binned[xbin]:
            occur[xbin]+=1
            if compare_to:
                compval=compare_to(k)
                theor[xbin]+= compval
                comp.setdefault( (compval,),0.)
                comp[ (compval,) ]+=dic[k]
    theor[occur!=0]/=occur[occur!=0]
    if compare_to:
        comp=binsumbase(comp,nbins=50)
    else:
        comp=None,None
    return bins,binval,theor, comp[0],comp[1]



def bayes(dic, col=0, nbins=20,stateproba=None):
    #Let x be one bin-value
    #P(Event|x)=P(x|Event)*P(Event)/P(x)
    #P(x|Event)= n_keys(x)/n_keys_tot
    #P(Event)=sum_keys dic[key]
    #P(x)=  n_states(x)/n_states_tot
    
    p_event=sum(dic[k] for k in dic)/len(ni)/len(nj)
    
    bins,binned=binner(dic,col,nbins)
    result=[]
    totlike=sum(sum([dic[k] for k in binned[xbin]]) for xbin in binned)
    for x in range(len(bins)):
        xproba=bisct.bisect_left( stateproba[col][0],bins[x] )
        if xproba>len(stateproba[col][0])-1:
            print x, stateproba[col][0][-1]
            xproba=len(stateproba[col][0])-1
        p_x=stateproba[col][1][xproba]/sum(stateproba[col][1])
        likelihood=sum([dic[k] for k in binned.get(x,[])])/totlike
        res=likelihood*p_event/p_x
        result.append(res)

    return np.array(bins), np.array(result)

#================================================================================
#================================================================================
#================================================================================
#============================== ANALYTICAL DEFS =================================
#================================================================================
#================================================================================
#================================================================================

def found(key):
    const=constants[tuple(key[-2:])]
    f,s,b,no,ono,tauh,taul,d,i,j=key
    ts=max(1.,const["taustheo"]) #* (  const["GRD_mx2"]/sqrt(2.) )/ max(1.,d)  )
    #ts/=(1-float(b)/N)
    if USE_RANDOM:
        ts=max(1.,(1+np.random.exponential(1.)*max(0.,ts-1.)))
    searchers=f-s
    return searchers/ts
    
    taul,tauh,D=taucalc(const )
    searchradius=sqrt(D*ts)
    if d>0:
        searchers=max(1.,min(searchers,searchers*( d/searchradius )**2  ))
    return searchers/ts
    
def bound(key):
    ff=found(key)
    f=key[0]
    s=key[1]
    b=key[2]
    return (f-s-ff)*min(1.,ff*constants[tuple(key[-2:])]["PC_lambda"]**2)
    
def left(key,bgo):
    f,s,b,no,ono,tauh,taul,d,i,j=key
    if tauh >0 and taul >0:
        o=min(1.,bgo)
        return s*min(1.,max(1./taul,o/tauh))
    else:
        return 0
        

def get_td(key):
    f,s,b,no,ono,tauh,taul,d,i,j=key
    const=constants[(i,j)]
    N=const["PC_n"]
    td=d/const["PC_v"]
    if USE_RANDOM:
        td=max(1,1+(td-1)*(np.random.exponential(1.)))
    return td

def lost(key,bgo,d=1.):
    f,s,b,no,ono,tauh,taul,d,i,j=key
    const=constants[(i,j)]
    if s>0 and no>0:
        res=left(key,bgo)*b/s
        #Time to diffuse away
        #taul,tauh,D=taucalc(const )
        #tdif= d**2/2./ D 
        #res/=max(1.,tdif)
        return min(b,res)
    else:
        return 0

    
def taucalc(const):
    th= const["PF_n"]/const["PC_q"]
    PS_p=const["PS_p"]
    if PS_p>0:
        tl= 1./PS_p #average time between jumps
    else:
        tl= 10000000000
    D=const["PC_v"]**2/const["PC_rp"]
    return tl, th,D

def discretize(variables):
    res=[]
    for v in variables:
        newv=int(round(v))
        if np.random.random()<(v%1.0):
            newv+=1
        res.append(newv)
    return res
    
def dcalc(bs,D,key,const):
    #Distance between fishers
    
    f,s,b,no,ono,tauh,taul,d,i,j=key
    if f+b<=1:
        return 0
    return const["GRD_mx2"] /sqrt(2.)
    
    if 0:
        #Random variables method
        nedges=N*(N-1)/2
        
        d= hypot(*np.random.uniform(-1.,1.,2)*const["GRD_mx2"]) *(nedges-bs) #+ (d/2.)*(b-bs)
        #for i in range(int(f-s)):
        #for i in range(int(b-bs)):
        return d/nedges
    
    #ODE Method
    diff=max(.01,(f-s)/const["PC_n"]*D/max(1.,d ))
    
    contract=max(.01,(b-bs)/const["PC_n"])#*const["PC_v"])
    
    dd=np.random.normal(diff/2.,.5)  -np.random.normal(contract,.5)
    d+=dd
    if d>const["GRD_mx2"]/sqrt(2.):
        d=2*const["GRD_mx2"]/sqrt(2.)-d
    if d<0:
        d=0
    return d 
    
def simu(fsource,length=100000,use_random=USE_RANDOM):
    i,j=fsource
    const=constants[(i,j)]
    t=0
    N=const["PC_n"]
    Sn=const["PS_n"]
    ncliq=const["PC_ncliq"]
    if RANDOM_PARTITION:
        cliquesizes=random_partition(N)
        f=np.array(cliquesizes)
        #print 'cliques:',f
        ncliq=len(f)
    else:    
        f=np.ones(ncliq)*floor(N/ncliq) #Free fishers
    s=np.zeros(ncliq) #Free fishers in a school
    b=np.zeros(ncliq) #Bound fishers
    bs=np.zeros(ncliq) #Bound fishers in a school
    no=np.zeros(ncliq) #Number of occupied schools
    ono=np.zeros(ncliq) #Total number of fishers in a school (<o>*n_o)
    taul,tauh,D= taucalc(const )
    res=[]
    
    d=np.ones(ncliq)* const["GRD_mx2"]/2.

    key=np.array((f,s,b,bs,no, np.ones(ncliq)* tauh, np.ones(ncliq)* taul, np.ones(ncliq)* d, np.ones(ncliq)* i, np.ones(ncliq)* j))

    dd=[]

    events={}
    def addevent(tup):
        e,ekey=tup
        ekey=tuple(ekey)
        events.setdefault(e,{})
        events[e].setdefault(ekey,0.)
        events[e][ekey]+=1.

    tsteps=4
    deltat=1./tsteps

    step=0
    bysize=np.zeros( (N,2))
    
    while t<length:
        step+=1
        t+=deltat
        #bgo=np.sum(key[1]+key[3])/np.sum(key[4])# BGO Meanfield version
        for cliq in range(ncliq):
            #BGO: not completely mean-field version (i.e. sum separately on this clique and others)
            tmpkey=np.ma.array(key, mask=False)
            tmpkey.mask[:,cliq]=True
            if np.sum(tmpkey[4])>0:
                bgo=np.sum(tmpkey[1]+tmpkey[3])/np.sum(tmpkey[4])
            else:
                bgo=0
            if key[4,cliq]>0:
                bgo+=(key[1,cliq]+key[3,cliq])/key[4,cliq]
        
            f,s,b,bs,no,th,tl,d,i,j=key[:,cliq]
            taus=constants[tuple(key[:,cliq][-2:])]["taustheo"]
            Wfs=found(key[:,cliq])
            Wfb=bound(key[:,cliq])
            Wbf=lost(key[:,cliq],bgo)
            Wsf=left(key[:,cliq],bgo)


            Wfs*=deltat
            Wfb*=deltat
            Wsf*=deltat
            Wbf*=deltat

            if DISCRETIZE_SIMU:
                Wfs,Wfb,Wbf,Wsf=discretize([Wfs,Wfb,Wbf,Wsf])
                Wfs=min(f-s,Wfs)
                Wfb=min(max(0,f-s-Wfs-1),Wfb)
                Wbf=min(b,Wbf)
                Wsf=min(s,Wsf)
                for x in range(int(Wfs)):
                    addevent(("simu_found_school",key[:,cliq]))
                for x in range(int(Wfb)):
                    addevent(("simu_bound",key[:,cliq]))
                for x in range(int(Wbf)):
                    addevent(("simu_lost_school",key[:,cliq]))
                for x in range(int(Wsf)):
                    addevent(("simu_left_school",key[:,cliq]))
            else:
                Wfs=min(f-s,Wfs)
                Wfb=min(max(0,f-s-1),Wfb)
                Wbf=min(b,Wbf)
                Wsf=min(s,Wsf)

            assert s<=f
            assert not np.isnan(s)
            df=Wbf-Wfb#+Wsf-Wfs-Wfb
            ds=Wfs-Wsf
            db=Wfb-Wbf
            #dno=Wfs-Wsf
            
            f+=df
            s+=ds
            b+=db
            #no+=dno
            no=Sn*(1.-exp(-s/Sn))
            d=dcalc(bs,D,key[:,cliq],const)
            dd.append(d)

            if DISCRETIZE_SIMU:       
                no= int(round(no))
            #Typical distance between fishers
            td=get_td(key[:,cliq])
            if b>0 and s>0:
                dbs=(b-bs)/td - Wsf*bs/s
                bs+=dbs*deltat
                if bs<0:
                    bs=0
            else:
                bs=0
            ono=s+bs

            if use_random:
                th,tl=tauh,1+(taul-1)*np.random.exponential(1.) #1+(tauh-1)*np.random.exponential(1.)
            else:
                th,tl=tauh,taul
            key[:,cliq]=(f,s,b,bs,no,th,tl,d,i,j)
            if RANDOM_PARTITION:
                bysize[cliquesizes[cliq]-1]+= ((float(bs)+s)/cliquesizes[cliq],1)

        
        if step%tsteps==0:
            res.append(np.sum(key,axis=1) )
    if RANDOM_PARTITION:
        bysize[bysize==0]=np.NaN
        bysize=bysize[:,0]/bysize[:,1]
    else:
        bysize=bysize[:,0]
        
        
            
    if 0:
        hist(np.array(dd)/const["GRD_mx2"]*sqrt(2.) ,hold=1,normed=1,log='y',bins=200)    
        plot(*distdic[(i,j)])
    return np.array(res),events,bysize

#================================================================================
#================================================================================
#================================================================================
#========================   DO SIMULATION  =============================
#================================================================================
#================================================================================
#================================================================================



H=np.zeros((len(ni),len(nj),2))
H_by_size=np.zeros((len(ni),len(nj),max( constants[k]["PC_n"] for k in constants  )) )
math=np.zeros((len(ni),len(nj),4))
mathsimu=np.zeros((len(ni),len(nj),4))
simustates={}
ix=-1
for i in ni:
    print i
    ix+=1
    iy=-1
    for j in nj:
        const=constants[(i,j)]
        iy+=1
        s=series[(i,j)]
        simulated,simuevents,bysize=[],[],[]
        for ens in range(nens):
            sim,simevt,bsize =simu((i,j),s.shape[0]/nens*3)
            simulated.append(sim)
            simuevents.append(simevt)
            bysize.append(bsize)
        print '\n'
        serieslength=min(len(ser) for ser in simulated)
       # simulated=np.mean([ser[:serieslength] for ser in simulated],0)
        simulated=np.concatenate(simulated)
        bysize = np.ma.masked_array(bysize,np.isnan(bysize))
        bysize=np.mean(bysize,0).filled(np.nan)*const["PC_q"]
        if RANDOM_PARTITION:
            print bysize        
            plot(cliqH[(i,j)],bysize,xs=range(len(bysize)))    
        H_by_size[ix,iy][:len(bysize)]=bysize
        Hrate=np.mean(s[:,1]+s[:,3])/const["PC_n"]*const["PC_q"]
        Hrate_simu=np.mean(simulated[:,1]+simulated[:,3])/const["PC_n"]*const["PC_q"]
        H[ix,iy]=Hrate,Hrate_simu
        math[ix,iy]=np.mean(s[:,1]),np.mean(s[:,2]),np.mean(s[:,3]),np.mean(s[:,4])
        mathsimu[ix,iy]=np.mean(simulated[:,1]),np.mean(simulated[:,2]),np.mean(simulated[:,3]),np.mean(simulated[:,4])
        

        for x in range(simulated.shape[0]):
            key=tuple(round(x) for x in simulated[x] )
            simustates.setdefault(key,0)
            simustates[key]+=1.        
            
        for simevt in simuevents:
            for e,vals in simevt.iteritems() :
                evtdic.setdefault(e,{})
                for key in vals:    
                    evtdic[e].setdefault(key,0.)
                    evtdic[e][key]+=vals[key]/s.shape[0]

        
        if 0:
            plt.subplot("211")
            hist(simulated[:,1],log='y',hold=1)
            hist(s[:,1],log='y',hold=1)
            plt.subplot("212")
            hist(simulated[:,2],log='y',hold=1)
            hist(s[:,2],log='y')
        if 0:
            plot( runavg(s[:,1],window)/  runavg(s[:,0],window),title='Bound/Finder')
        window=5#50

        if 0:
            #RUNAVG
            plt.subplot("221")
            print np.mean(simulated[:,3]),np.mean(s[:,3]), s.shape
            plot( runavg(simulated[:,1],window) ,hold=1)
            plot( runavg(s[:,1],window),title="s" ,hold=1,alpha=.5)
            plt.subplot("222")
            plot( runavg(simulated[:,2],window) ,hold=1)
            plot( runavg(s[:,2],window),title="b" ,hold=1,alpha=.5)
            plt.subplot("223")
            plot( runavg(simulated[:,3],window) ,hold=1)
            plot( runavg(s[:,3],window),title="bs" ,hold=1,alpha=.5)
            plt.subplot("224")
            plot( runavg(simulated[:,4],window) ,hold=1)
            plot( runavg(s[:,4],window),title="no",alpha=.5)
        if 0:
            #CUMAVG
            plt.subplot("221")
            plot( cumavg(simulated[:,1]) ,hold=1)
            plot( cumavg(s[:,1]),title="s" ,hold=1,alpha=.5)
            plt.subplot("222")
            plot( cumavg(simulated[:,2]) ,hold=1)
            plot( cumavg(s[:,2]),title="b" ,hold=1,alpha=.5)
            plt.subplot("223")
            plot( cumavg(simulated[:,3]) ,hold=1)
            plot( cumavg(s[:,3]),title="bs" ,hold=1,alpha=.5)
            plt.subplot("224")
            plot( cumavg(simulated[:,4]) ,hold=1)
            plot( cumavg(s[:,4]),title="no",alpha=.5)




print "Making real stateproba"
realstateproba=[binsumbase(states,col) for col in (1,2)]

print "Making simulation stateproba"
simustateproba=[binsumbase(simustates,col) for col in (1,2)]

if 1:
    plot(* realstateproba[0],hold=1)
    plot(* simustateproba[0],color='g',title='Histogram: finders')

    plot(* realstateproba[1],hold=1)
    plot(* simustateproba[1],color='g',title='Histogram: bound')

#================================================================================
#================================================================================
#================================================================================
#=============================   FINAL PLOTS   ==================================
#================================================================================
#================================================================================
#================================================================================


if 1:
    #Catch versus theory
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_zlabel(r"$H$")

    xkey=sorted(x_ax)[0]
    ykey="PS_n"#sorted(y_ax)[0]
    X=x_ax[xkey]
    Y=y_ax[ykey]
    plt.xlabel(xkey)
    plt.ylabel(ykey)
    
    ax.scatter(X.ravel(),Y.ravel(),H[:,:,0].ravel())
    ax.scatter(X.ravel(),Y.ravel(),data["Hrate"].ravel(),color='r')
    ax.scatter(X.ravel(),Y.ravel(),H[:,:,1].ravel(),color='g')
    ax.plot_wireframe(X,Y,H[:,:,1],color='g')
    plt.show()
    for i in range(4):
        title=('finders','bound','bs','no')
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(X.ravel(),Y.ravel(),math[:,:,i].ravel())
        ax.scatter(X.ravel(),Y.ravel(),mathsimu[:,:,i].ravel(),color='g')
        plt.title(title[i])    
        plt.show()    
    
def bayesplot(col=0,simu=False):
    pairs=[("left_school",left, "Left school"),("lost_school",lost, "Lost school"),("found_school",found, "Found school"),("bound",bound, "Heed call") ]
    coltitle=("Finders","Bound")
    if simu:
        stateproba=simustateproba
    else:
        stateproba=realstateproba
    for k,f,title in pairs:
        if simu:
            title="(SIMU) "+title
            k='simu_'+k
        bins,vals, theo,compbins,compvals=binsum(evtdic[k],col,compare_to= f,stateproba=stateproba)
        bins,vals=bayes(evtdic[k],col,stateproba=stateproba)
        plot(vals/np.mean(vals),theo/np.mean(theo),xs=bins,title=title,xlabel=coltitle[col])
    

if 1:
    print "BAYES FOR ABM (FINDERS)"
    bayesplot(1,0)

if 1:
    print "BAYES FOR SIMU (FINDERS) -- REQUIRES DISCRETIZE_SIMU"
    bayesplot(1,1)

if 1:    
    print "BAYES FOR ABM (BOUND)"
    bayesplot(2,0)

if 1:
    print "BAYES FOR SIMU (BOUND) -- REQUIRES DISCRETIZE_SIMU"
    bayesplot(2,1)
    
if 0:
    print "SCATTER COMP"

    bins,vals, theo,compbins,compvals=binsum(evtdic["found_school"],0,compare_to= found,stateproba=realstateproba)
    plot(vals,theo,xs=bins,title='Found school/finders')
    scatter(compbins,compvals,log='xy')

    bins,vals, theo,compbins,compvals=binsum(evtdic["bound"],0,compare_to= bound,stateproba=realstateproba)
    plot(vals,theo,xs=bins,title='Heed call/finders')
    scatter(compbins,compvals,log='xy')
    
    bins,vals, theo,compbins,compvals=binsum(evtdic["found_school"],1,compare_to= found,stateproba=realstateproba)
    scatter(compbins,compvals,log='xy')

    bins,vals, theo,compbins,compvals=binsum(evtdic["bound"],1,compare_to= bound,stateproba=realstateproba)
    scatter(compbins,compvals,log='xy')

