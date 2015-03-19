from datatools import *
from plot_benchcore import *

ni,nj=range(1,5),range(1,5)

USE_RANDOM=1

filename="Data_anaflux"
set_constants(filename)
data=load_data(filename) 

keynames=("s","b","no","o * no","tauh" ,"taul")

evtdic={e:{} for  e in ("found_school","left_school","lost_school", "bound")}

constants={}
states={} #States throughout the simulation
series={} #Time series of te states



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
        d=np.load("../Data/dist-{}-{}.npy".format(i,j))
        s=np.load("../Data/states-{}-{}.npy".format(i,j)) 
        
        series[(i,j)]=s
        
        for x in range(s.shape[0]):
            key=tuple(s[x] )
            key+=(i,j)
            states.setdefault(key,0)
            states[key]+=1.
            
        lines=open("../Data/turn{}-{}.dat".format(i,j),'r').readlines()
        nturns,ts,tstheo=[float(l.strip()) for l in lines[:3] ]
        constants[(i,j)]={l.split()[0]:eval(l.split()[1])  for l in lines[3:] }
        constants[(i,j)]["taus"]=ts
        constants[(i,j)]["taustheo"]=tstheo 
                   
        for  e in ("found_school","left_school","lost_school", "bound"):
            localdic={}
            arr=np.load("../Data/{}-{}-{}.npy".format(e,i,j))

            for x in range(arr.shape[0]):
                key=tuple(arr[x])
                key+=(i,j )
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

stateproba=[binsumbase(states,col) for col in (0,1)]

def binsum(dic, col=0,nbins=100,compare_to=None):
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



def bayes(dic, col=0, nbins=20):
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


def found(key):
    N=constants[key[-2:]]["PC_n"]
    s=key[0]
    b=key[1]
    f=N-b
    ts=constants[key[-2:]]["taustheo"]
    #ts/=(1-float(b)/N)
    if  USE_RANDOM:
        ts=(1+np.random.exponential(1.)*max(0.,ts-1.))
    return (f-s)/ts
    
def bound(key):
    ff=found(key)
    N=constants[key[-2:]]["PC_n"]
    s=key[0]
    b=key[1]
    f=N-b
    return (f-s-ff)*min(1.,ff*constants[key[-2:]]["PC_lambda"]**2)
    
def left(key):
    s,b,no,ono,tauh,taul,i,j=key
    if tauh >0 and taul >0 and no>0:
        o=1.#min(1.,ono/no)
        return s*max(1./taul,o/tauh)
    else:
        return 0
        

def lost(key):
    s,b,no,ono,tauh,taul,i,j=key
    if s>0 and no>0:
        return left(key)*b/s
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
    
def simu(fsource,length=100000,use_random=USE_RANDOM):
    i,j=fsource
    const=constants[(i,j)]
    t=0
    N=const["PC_n"]
    f=N
    s=0
    b=0
    no=0 #Number of occupied schools
    ono=0 #Total number of fishers in a school (<o>*n_o)
    bs=0 #Bound fishers in a school
    taul,tauh,D=taucalc(const )
    res=[]
    key=(s,b,no,ono,tauh,taul,i,j)
    while t<length:
        t+=1
        taus=constants[key[-2:]]["taustheo"]
        Wfs=found(key)
        Wfb=bound(key)
        Wbf=lost(key)
        Wsf=left(key)
        f+=Wbf+Wsf-Wfs-Wfb
        s+=Wfs-Wsf
        b+=Wfb-Wbf
        no+=Wfs-Wsf
        #Typical distance between fishers
        if no>0 and ono>0:
            d=min(const["GRD_mx2"]/sqrt(2.),sqrt(D/no)*N/(ono/no))
        else:
            d=const["GRD_mx2"]/sqrt(2.)
        td=d/const["PC_v"]
        if use_random:
            td=max(1,td*(.5+np.random.exponential(.5)))
        if b>0:
            bs+=(b-bs)/td - Wbf*bs/b
            if bs<0:
                bs=0
        else:
            bs=0
        ono=s+bs

        if use_random:
            th,tl=tauh,1+(taul-1)*np.random.exponential(1.) #1+(tauh-1)*np.random.exponential(1.)
        else:
            th,tl=tauh,taul
        key=(s,b,no,ono,th,tl,i,j)
        res.append(key)
        if np.isnan(s):
            break
    return np.array(res)

H=np.zeros((len(ni),len(nj),2))
math=np.zeros((len(ni),len(nj),4))
mathsimu=np.zeros((len(ni),len(nj),4))
ix=-1
for i in ni:
    ix+=1
    iy=-1
    for j in nj:
        iy+=1
        s=series[(i,j)]
        simulated=simu((i,j),s.shape[0])
        Hrate=np.mean(s[:,3])/constants[(i,j)]["PC_n"]*constants[(i,j)]["PC_q"]
        Hrate_simu=np.mean(simulated[:,3])/constants[(i,j)]["PC_n"]*constants[(i,j)]["PC_q"]
        H[ix,iy]=Hrate,Hrate_simu
        math[ix,iy]=np.mean(s[:,0]),np.mean(s[:,1]),np.mean(s[:,2]),np.mean(s[:,3])
        mathsimu[ix,iy]=np.mean(simulated[:,0]),np.mean(simulated[:,1]),np.mean(simulated[:,2]),np.mean(simulated[:,3])
        if 0:
            plt.subplot("211")
            hist(simulated[:,0],log='y',hold=1)
            hist(s[:,0],log='y',hold=1)
            plt.subplot("212")
            hist(simulated[:,1],log='y',hold=1)
            hist(s[:,1],log='y')
        window=50#50

        if 1:
            #RUNAVG
            plt.subplot("221")
            print np.mean(simulated[:,3]),np.mean(s[:,3])
            plot( runavg(simulated[:,0],window) ,hold=1)
            plot( runavg(s[:,0],window),title="s" ,hold=1,alpha=.5)
            plt.subplot("222")
            plot( runavg(simulated[:,1],window) ,hold=1)
            plot( runavg(s[:,1],window),title="b" ,hold=1,alpha=.5)
            plt.subplot("223")
            plot( runavg(simulated[:,2],window) ,hold=1)
            plot( runavg(s[:,2],window),title="no" ,hold=1,alpha=.5)
            plt.subplot("224")
            plot( runavg(simulated[:,3],window) ,hold=1)
            plot( runavg(s[:,3],window),title="ono",alpha=.5)
        if 0:
            #CUMAVG
            plt.subplot("221")
            plot( cumavg(simulated[:,0]) ,hold=1)
            plot( cumavg(s[:,0]),title="s" ,hold=1,alpha=.5)
            plt.subplot("222")
            plot( cumavg(simulated[:,1]) ,hold=1)
            plot( cumavg(s[:,1]),title="b" ,hold=1,alpha=.5)
            plt.subplot("223")
            plot( cumavg(simulated[:,2]) ,hold=1)
            plot( cumavg(s[:,2]),title="no" ,hold=1,alpha=.5)
            plt.subplot("224")
            plot( cumavg(simulated[:,3]) ,hold=1)
            plot( cumavg(s[:,3]),title="ono",alpha=.5)

        if 0:
            plot( runavg(s[:,1],window)/  runavg(s[:,0],window),title='Bound/Finder')

        for  e in ("simu_found_school","simu_left_school","simu_lost_school", "simu_bound"):
            arr=simulated
            evtdic.setdefault(e,{})

            for x in range(arr.shape[0]):
                key=tuple(arr[x])
                evtdic[e].setdefault(key,0.)
                evtdic[e][key]+=1./nturns

if 1:
    #Catch versus theory
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_zlabel(r"$H$")

    xkey=sorted(x_ax)[0]
    ykey=sorted(y_ax)[0]
    X=x_ax[xkey]
    Y=y_ax[ykey]
    plt.xlabel(xkey)
    plt.ylabel(ykey)
    
    ax.scatter(X.ravel(),Y.ravel(),H[:,:,0].ravel())
    ax.scatter(X.ravel(),Y.ravel(),data["Hrate"].ravel(),color='r')
    ax.scatter(X.ravel(),Y.ravel(),H[:,:,1].ravel(),color='g')
    ax.plot_wireframe(X,Y,H[:,:,1],color='g')
    plt.show()
    #Finders
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X.ravel(),Y.ravel(),math[:,:,0].ravel())
    ax.scatter(X.ravel(),Y.ravel(),mathsimu[:,:,0].ravel(),color='g')
    plt.show()
    #Bound
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X.ravel(),Y.ravel(),math[:,:,1].ravel())
    ax.scatter(X.ravel(),Y.ravel(),mathsimu[:,:,1].ravel(),color='g')
    plt.show()
    #no
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X.ravel(),Y.ravel(),math[:,:,2].ravel())
    ax.scatter(X.ravel(),Y.ravel(),mathsimu[:,:,2].ravel(),color='g')
    plt.show()

if 0:
    plot(stateproba[0][0],stateproba[0][1]/sum(stateproba[0][1])*6,hold=1)
    hist([k[0] for k in evtdic["found_school"]],normed=1,hold=1)
    hist([k[0] for k in evtdic["left_school"]],normed=1,hold=1)
    hist([k[0] for k in evtdic["bound"]],normed=1)
    plot(stateproba[1][0],stateproba[1][1]/sum(stateproba[1][1])*10,hold=1)
    hist([k[1] for k in evtdic["found_school"]],normed=1,hold=1)
    hist([k[1] for k in evtdic["left_school"]],normed=1,hold=1)
    hist([k[1] for k in evtdic["bound"]],normed=1)
    
    
if 0:
    bins,vals, theo,compbins,compvals=binsum(evtdic["found_school"],0,compare_to= found)
    plot(vals,theo,xs=bins,title='Found school/finders')
    scatter(compbins,compvals,log='xy')


    bins,vals, theo,compbins,compvals=binsum(evtdic["bound"],0,compare_to= bound)
    plot(vals,theo,xs=bins,title='Heed call/finders')
    scatter(compbins,compvals,log='xy')
    
    
if 0:
    
    bins,vals, theo,compbins,compvals=binsum(evtdic["simu_left_school"],0,compare_to= left)
    bins,vals=bayes(evtdic["simu_left_school"],0)
    plot(vals,theo,xs=bins,title='Bayes: Left school/finders')

    bins,vals, theo,compbins,compvals=binsum(evtdic["simu_lost_school"],0,compare_to= lost)
    bins,vals=bayes(evtdic["simu_lost_school"],0)
    plot(vals,theo,xs=bins,title='Bayes: Lost school/finders')
    
    bins,vals, theo,compbins,compvals=binsum(evtdic["simu_found_school"],0,compare_to= found)
    bins,vals=bayes(evtdic["simu_found_school"],0)
    plot(vals,theo,xs=bins,title='Bayes: Found school/finders')

    bins,vals, theo,compbins,compvals=binsum(evtdic["simu_bound"],0,compare_to= bound)
    bins,vals=bayes(evtdic["simu_bound"],0)
    plot(vals,theo,xs=bins,title='Bayes: Heed call/finders')
if 1:
    
    bins,vals, theo,compbins,compvals=binsum(evtdic["left_school"],0,compare_to= left)
    bins,vals=bayes(evtdic["left_school"],0)
    plot(vals,theo,xs=bins,title='Bayes: Left school/finders')

    bins,vals, theo,compbins,compvals=binsum(evtdic["lost_school"],0,compare_to= lost)
    bins,vals=bayes(evtdic["lost_school"],0)
    plot(vals,theo,xs=bins,title='Bayes: Lost school/finders')
    
    bins,vals, theo,compbins,compvals=binsum(evtdic["found_school"],0,compare_to= found)
    bins,vals=bayes(evtdic["found_school"],0)
    plot(vals,theo,xs=bins,title='Bayes: Found school/finders')

    bins,vals, theo,compbins,compvals=binsum(evtdic["bound"],0,compare_to= bound)
    bins,vals=bayes(evtdic["bound"],0)
    plot(vals,theo,xs=bins,title='Bayes: Heed call/finders')
    

if 1:    
    bins,vals, theo,compbins,compvals=binsum(evtdic["left_school"],1,compare_to= left)
    bins,vals=bayes(evtdic["left_school"],1)
    plot(vals,theo,xs=bins,title='Bayes: Left school/bound')

    bins,vals, theo,compbins,compvals=binsum(evtdic["lost_school"],1,compare_to= lost)
    bins,vals=bayes(evtdic["lost_school"],1)
    plot(vals,theo,xs=bins,title='Bayes: Lost school/bound')


    bins,vals, theo,compbins,compvals=binsum(evtdic["found_school"],1,compare_to= found)
    bins,vals=bayes(evtdic["found_school"],1)
    plot(vals,theo,xs=bins,title='Bayes: Found school/bound')

    bins,vals, theo,compbins,compvals=binsum(evtdic["bound"],1,compare_to= bound)
    bins,vals=bayes(evtdic["bound"],1)
    plot(vals,theo,xs=bins,title='Bayes: Heed call/bound')

if 0:

    bins,vals, theo,compbins,compvals=binsum(evtdic["found_school"],1,compare_to= found)
    scatter(compbins,compvals,log='xy')

    bins,vals, theo,compbins,compvals=binsum(evtdic["bound"],1,compare_to= bound)
    scatter(compbins,compvals,log='xy')

