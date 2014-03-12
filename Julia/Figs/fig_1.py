import numpy as np
import matplotlib.pyplot as plt

## Load in data
Fish = np.load("../Data/Data_fish.npy")
Cons = np.load("../Data/Data_fishers.npy")
Clus = np.load("../Data/Data_clusters.npy")

## plot
x   = Fish[:,0,100];
y   = Fish[:,1,100];
xc  = Cons[:,0,100];
yc  = Cons[:,1,100];
xl  = Clus[:,0,100];
yl  = Clus[:,1,100];

fig = plt.figure(1, figsize=(10, 8),edgecolor=[.4,.4,.4]);
plt.plot(x,y,'o',alpha=0.5,mec='none',color=[.3,.3,1.],ms=9);
plt.plot(xc,yc,'o',alpha=.75,mec='none',color=[1,0,0],ms=9);
plt.plot(xl,yl,'^',alpha=1,mec='none',color=[0,1,0],ms=9);

plt.plot(xc[1],yc[1],'ro',alpha=1,mec='k',mfc='none',ms=12);
plt.plot([xc[1],xc[2]],[yc[1],yc[2]],'-r',lw=2);
plt.plot([xc[1],xc[3]],[yc[1],yc[3]],'-r',lw=1);
plt.plot([xc[1],xc[4]],[yc[1],yc[4]],'-r',lw=.5);

plt.axis([0,100,0,100])
plt.xticks([]);
plt.yticks([]);

# save
plt.savefig('./PNG/Fig_1.png',dpi=600,bbox_inches='tight')
plt.close()


