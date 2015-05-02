from datatools import *
from plot_benchcore import *


filename="Data_LHFpopIFQ"
set_constants(filename)
data=load_data(filename) 

varx,vary='pop','quota'

xs=data[varx][0]
ys=data[vary].T[0]
H=data['Hrate']

try:
    LAM=data["PC_lambda"]
except:
    LAM=np.ones(H.shape)*constants["PC_lambda"]
FISH=data['pop']

start=(xs[-4],ys[-4])

def interpolate(pos,surf):
    ix=bisct.bisect_left(list(xs),pos[0])
    iy=bisct.bisect_left(list(ys),pos[1])
    ds=[ (0,0), (0,1), (1,0), (1,1) ]
    return sum(array([surf[ix+dx,iy+dy] for dx,dy in ds])/len(ds) )
    
def maximize(pos,surf,axis=1):
    ix=bisct.bisect_left(list(xs),pos[0])
    iy=bisct.bisect_left(list(ys),pos[1])
    if axis==1:
        new =np.argmax(surf[ix])
    if axis==0:
        new=np.argmax(surf[:,iy])
    return new
    
def bisect(traj):
    return [(bisct.bisect_left(list(xs),pos[0]),bisct.bisect_left(list(ys),pos[1])) for pos in traj]

def convert_pos(pos):
    newpos=pos.copy()
    if varx=='PF_n':
        newpos[0]/=constants["PS_n"]
    if vary!='PC_lambda':
        #do something
        pass
    return newpos
    
def traj(init,**kwargs):
    adapt_lambda=kwargs.get('adapt_lambda',False)
    
    fish=interpolate(init,FISH)
    lam=interpolate(init,LAM)
    t=0
    pos=array((fish,lam))
    poses=[]
    while t < 10:
        epos=convert_pos(pos)
        pos[0]-= interpolate(epos,H)*constants['PC_n']
        if adapt_lambda:
            pos[1]=maximize(epos,H,1)
        t+=1
        poses.append(pos.copy())
    return array(poses)

t1=traj(start).T
#plot(*t1,hold=1)
print bisect(t1)
scatter(t1[0],t1[1],H[bisect(t1)])
