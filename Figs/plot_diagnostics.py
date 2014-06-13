import numpy as np
import matplotlib.pyplot as plt

## Load in data
SN      = np.load("../Data/Data_evo_SN.npy")
cpue_mu = np.load("../Data/Data_evo_cpue_mu.npy")
cpue_s2 = np.load("../Data/Data_evo_cpue_s2.npy")

## plot
x  = np.arange(0,SN.shape[2]);
y1 = np.mean(SN,1);
y2 = np.mean(y1,0);

y3 = cpue_mu;
y4 = np.mean(y3,0);

y5 = cpue_s2;
y6 = np.mean(y5,0);

## Plot Connectance over evo time
fig, ax = plt.subplots(1, figsize=(10, 8),edgecolor=[.4,.4,.4]);
for i in np.arange(0,y1.shape[0]):
    plt.plot(x,y1[i,:],'-',alpha=.5,color=[0,0,1],lw=1);
plt.plot(x,y2,'-',alpha=1,color=[1,0,0],lw=2);
ax.set_ylabel('Connectance', color='k')
ax.set_xlabel('Time', color='k')
#plt.savefig('./PNG/Fig_diag_1.png',dpi=600,bbox_inches='tight')


## Plot mean CPUE over evo time
fig, ax = plt.subplots(1, figsize=(10, 8),edgecolor=[.4,.4,.4]);
for i in np.arange(0,y1.shape[0]):
    plt.plot(x,y3[i,:],'-',alpha=.5,color=[0,0,1],lw=1);
plt.plot(x,y4,'-',alpha=1,color=[1,0,0],lw=2);
ax.set_ylabel('Expected catch per unit time', color='k')
ax.set_xlabel('Time', color='k')
#plt.savefig('./PNG/Fig_diag_2.png',dpi=600,bbox_inches='tight')


## Plot variance in CPUE over evo time
fig, ax = plt.subplots(1, figsize=(10, 8),edgecolor=[.4,.4,.4]);
for i in np.arange(0,y1.shape[0]):
    plt.plot(x,y5[i,:],'-',alpha=.5,color=[0,0,1],lw=1);
plt.plot(x,y6,'-',alpha=1,color=[1,0,0],lw=2);
ax.set_ylabel('Variance in catch per unit time', color='k')
ax.set_xlabel('Time', color='k')
#plt.savefig('./PNG/Fig_diag_3.png',dpi=600,bbox_inches='tight')



plt.show()
