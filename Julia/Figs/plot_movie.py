import numpy as np
import matplotlib.pyplot as plt

## Load in data
Fish   = np.load("../Data/Data_fish.npy")
Cons = np.load("../Data/Data_fishers.npy")

## plot

for i in np.arange(0,Fish.shape[2]):

    x   = Fish[:,0,i];
    y   = Fish[:,1,i];
    xc  = Cons[:,0,i];
    yc  = Cons[:,1,i];

    fig = plt.figure(1, figsize=(10, 8),edgecolor=[.4,.4,.4]);
    plt.plot(x,y,'o',alpha=0.5,mec='none',color=[.3,.3,1.],ms=9);
    plt.plot(xc,yc,'o',alpha=1,mec='none',color=[1,0,0],ms=9);
    plt.axis([0,100,0,100])
    plt.xticks([]);
    plt.yticks([]);

    # save
    ti = str(i+100000)
    plt.savefig('./PNG/Fig_' +ti[1:] +'.png',dpi=100,bbox_inches='tight')
    print(i)
    plt.close('all')


