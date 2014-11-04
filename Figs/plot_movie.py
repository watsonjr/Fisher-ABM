import numpy as np
import matplotlib.pyplot as plt
from subprocess import call
import os



imgformat='jpg'
maxpic=500 #Maximum number of pictures made


## Remove old figs
try:
    call("rm -f ./Movie/*.{}".format(imgformat),shell=True )
except:
    print("Your OS is not supported by my awesome movie plotter.")

## Load in data
Fish = np.load("../Data/Data_fish.npy")
Cons = np.load("../Data/Data_fishers.npy")
Clus = np.load("../Data/Data_clusters.npy")
Pop = np.load("../Data/Data_cluspop.npy")
Harv = np.load("../Data/Data_harvest.npy")
MI = np.load("../Data/Data_MI.npy")

implicit_fish=True

#GET PARAMETER VALUES
fcst=open("../Data/Data_figs.dat","r")
for l in fcst:
    line=l.split()
    if len(line)>1:
        exec("{} = {}".format(line[0],line[1]))
fcst.close()

#Size conversion (data to figure size)
def convert_dist(ax,d):
    return float(d)/GRD_mx* ( ax.transData.transform((0,GRD_mx ) )-ax.transData.transform((0,0) ) )[1]

## plot

for i in np.arange(0,min(maxpic,Fish.shape[2])):

    x   = Fish[:,0,i];
    y   = Fish[:,1,i];
    xc  = Cons[:,0,i];
    yc  = Cons[:,1,i];
    xl  = Clus[:,0,i];
    yl  = Clus[:,1,i];

    fig = plt.figure( edgecolor=[.4,.4,.4]);
    ax = fig.add_subplot(111)
    ax.set(aspect=1)
    plt.axis([0,GRD_mx,0,GRD_mx])
    
    if not implicit_fish:
        plt.plot(x,y,'o',alpha=0.5,mec='none',color=[.3,.3,1.],ms=9);
        plt.plot(xl,yl,'^',alpha=1,mec='none',color=[0,1,0],ms=9);
    else:
        plt.scatter(xl,yl,s=3.14*convert_dist(ax,PF_sig)**2* Pop[:,i]/PF_n,alpha=.4 );
    plt.scatter(xc,yc,alpha=.3,color=[1,0,0],s=3.14*convert_dist(ax,PC_f)**2);
    #plt.scatter(xc,yc,alpha=1,c=[[1,(ca%PF_n)/PF_n,(ca%PF_n)/PF_n] for ca in Harv[:,i]],s=3.14*convert_dist(ax,PC_h)**2);
    plt.scatter(xc,yc,alpha=1,c=[[1,0,0] if mi else [1,1,1] for mi in MI[:,i]],s=3.14*convert_dist(ax,PC_h)**2);
    plt.xticks([]);
    plt.yticks([]);

    # save
    ti = str(i+100000)
    plt.savefig('./Movie/Fig_' +ti[1:] +'.{}'.format(imgformat),dpi=100,bbox_inches='tight')
    print(i)
    plt.close('all')

try:
    call("ffmpeg -y -r 18  -i Movie/Fig_%05d.{} movie.mp4".format(imgformat),shell=True)
except:
    pass
