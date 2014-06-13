import numpy as np
import matplotlib.pyplot as plt

## Load in data
MU = np.load("../Data/Data_opt_cpue_mu.npy")
S2 = np.load("../Data/Data_opt_cpue_s2.npy")
S = np.load("../Data/Data_opt_S.npy")

## plot
x   = S;
y1  = np.mean(MU,0);
y2  = np.mean(np.sqrt(S2),0);

fig, ax = plt.subplots(1, figsize=(10, 8),edgecolor=[.4,.4,.4]);
ax.plot(x,y1,'-',alpha=1,mec='none',color=[.3,.3,1.],lw=2);
ax.set_ylabel('Expected catch per season', color='k')
ax.set_xlabel('Pro-sociality', color='k')

fig, ax = plt.subplots(1, figsize=(10, 8),edgecolor=[.4,.4,.4]);
ax.plot(x,y2,'-',alpha=1,mec='none',color=[1,.3,.3],lw=2);
ax.set_ylabel('STD catch per season', color='k')
ax.set_xlabel('Pro-sociality', color='k')

plt.show()


# save
#plt.savefig('./PNG/Fig_1.png',dpi=600,bbox_inches='tight')


