from datatools import *
from plot_benchcore import *


constants={}


def found(key):
    N=constants[key[-2:]]["PC_n"]
    s=key[0]
    b=key[1]
    f=N-b
    ts=constants[key[-2:]]["taustheo"]
    return (f-s)/ts
    
def bound(key):
    ff=found(key)
    N=constants[key[-2:]]["PC_n"]
    s=key[0]
    b=key[1]
    f=N-b
    return (f-s-ff)*ff*constants[key[-2:]]["PC_lambda"]**2
    
def left(key):
    s,b,no,ono,tauh,taul,i,j=key
    if tauh >0 and taul >0 and no>0:
        o=min(1.,ono/no)
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
    
def simu(fsource,length=100000,use_random=1):
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
        if use_random:
            Wfs*=np.random.exponential(1.)
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
        td=d/const["PC_v"]*(1+np.random.exponential(1.))
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


dft={"GRD_mx2":50.,"PC_lambda":1.,"PC_v":1.,"PF_n":10,"PS_n":10,"PC_n":10,"PS_p":0.01, "PC_rp":.5, "PC_q":.1, "taustheo":100}

vari=("PC_lambda",np.linspace(.1,1.,4) )
varj=("PS_p",np.logspace(.005,.05,4))

ni,nj=range(0,4),range(0,4)
simulated={}
H=np.zeros((len(ni),len(nj)))
for i in nj:
    print i
    for j in nj:
        constants[(i,j)]={}
        constants[(i,j)].update(dft)
        constants[(i,j)][vari[0]]=vari[1][i]
        constants[(i,j)][varj[0]]=varj[1][j]
        simulated[(i,j)]=s=simu((i,j))
        Hrate=np.mean(s[:,3])/constants[(i,j)]["PC_n"]*constants[(i,j)]["PC_q"]
        H[i,j]=Hrate



#Catch
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.xlabel(r"X")
plt.ylabel(r"Y")
ax.set_zlabel(r"$H$")
X=np.array(list(ni)*len(nj) ).reshape((len(ni),len(nj)))
Y=np.array(list(nj)*len(ni) ).reshape((len(ni),len(nj))).T
ax.scatter(X.ravel(),Y.ravel(),H[:,:].ravel(),color='g')
ax.plot_wireframe(X,Y,H[:,:],color='g')
plt.show()

