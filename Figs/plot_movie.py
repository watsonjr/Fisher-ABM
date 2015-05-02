import numpy as np
import matplotlib.pyplot as plt
from subprocess import call
import os



imgformat='jpg'
maxpic=900 #Maximum number of pictures made
in_path="../Data"
#in_path="../LHFMovie2"

out_path="./Movie"
#out_path="./LHFMovie"

## Remove old figs
try:
    call("rm -f "+out_path+"/*.{}".format(imgformat),shell=True )
except:
    print("Your OS is not supported by my awesome movie plotter.")

## Load in data
Fish = np.load(in_path+"/Data_fish.npy")
Cons = np.load(in_path+"/Data_fishers.npy")
Clus = np.load(in_path+"/Data_clusters.npy")
Pop = np.load(in_path+"/Data_cluspop.npy")
Harv = np.load(in_path+"/Data_harvest.npy")
MI = np.load(in_path+"/Data_MI.npy")

implicit_fish=True

#GET PARAMETER VALUES
fcst=open(in_path+"/Data_figs.dat","r")
for l in fcst:
    line=l.split()
    if len(line)>1:
        exec("{} = {}".format(line[0],line[1]))
fcst.close()

#Size conversion (data to figure size)
def convert_dist(ax,d):
    return float(d)/GRD_mx* ( ax.transData.transform((0,GRD_mx ) )-ax.transData.transform((0,0) ) )[1]

## plot
maxHarv=np.amax(Harv)

for i in np.arange(0,min(maxpic,Fish.shape[2])):

    x   = Fish[:,0,i];
    y   = Fish[:,1,i];
    xc  = Cons[:,0,i];
    yc  = Cons[:,1,i];
    xl  = Clus[:,0,i];
    yl  = Clus[:,1,i];

    fig = plt.figure( edgecolor=[.4,.4,.4]);
    fig.set_size_inches(9.,6.)
    ax = plt.subplot2grid((2,2), (0, 0), rowspan=2)
    ax.set(aspect=1)
    plt.axis([0,GRD_mx,0,GRD_mx])
    
    if not implicit_fish:
        plt.plot(x,y,'o',alpha=0.5,mec='none',color=[.3,.3,1.],ms=9);
        plt.plot(xl,yl,'^',alpha=1,mec='none',color=[0,1,0],ms=9);
    else:
        #plt.scatter(xl,yl,s=3.14*convert_dist(ax,PF_sig)**2* Pop[:,i]/PF_n,alpha=.4 );
        ttt=np.array([.3,.3,1.])
        www=np.array([1.,1.,1.])
        #color=[www +(ttt-www) * Pop[j,0,i]/PF_n for j in range(PS_n)]
        plt.scatter(xl,yl,s=3.14*convert_dist(ax,PF_sig)**2,alpha=.4,color= [.3,.3,.5] );
        
    colors=[[1,0,0] if ca>=0 else [.8,.8,.8] for ca in MI[:,i]]
    plt.scatter(xc,yc,alpha=.3,color=colors,s=3.14*convert_dist(ax,PC_f)**2);
    
    #plt.scatter(xc,yc,alpha=1,c=[[1,(ca%PF_n)/PF_n,(ca%PF_n)/PF_n] for ca in Harv[:,i]],s=3.14*convert_dist(ax,PC_h)**2);

    #colors=[1,0,0]
    #colors=[[1,0,0] if mi else [1,1,1] for mi in MI[:,i]]
    #colors=[[1,1,1] if i>0 and Harv[f,i]!=Harv[f,i-1] else [1,0,0] for  f in range(PC_n)]
    
    plt.scatter(xc,yc,alpha=1,s=3.14*convert_dist(ax,PC_h)**2,c=colors);
    plt.xticks([]);
    plt.yticks([]);

    ax2 = plt.subplot2grid((2,2), (0,1))
    ax2.set_ylim(0,maxHarv )
    ax2.scatter(range(Harv.shape[0]),Harv[:,0,i])
    ax3 = plt.subplot2grid((2,2), (1, 1))
    ax3.plot([0]+[sum(Harv[:,0,j])-sum(Harv[:,0,j-1]) for j in range(1,i+1)])
    
    # save
    print(i)
    ti = str(i+100000)
    plt.savefig(out_path+'/Fig_' +ti[1:] +'.{}'.format(imgformat),dpi=100,bbox_inches='tight')
    plt.close('all')

try:
    call("ffmpeg -y -i "+out_path+"/Fig_%05d.{} -r 18  movie.mp4".format(imgformat),shell=True)
except:
    pass
