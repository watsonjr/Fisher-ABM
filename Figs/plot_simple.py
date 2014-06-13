import numpy as np
import matplotlib.pyplot as plt

## Load in data
Data = np.load("../Data/Data_simple.npz")
CPUE = Data['CPUE']
Tau  = Data['Tau']

## plot
x = np.linspace(0,1,10);
y1m = np.mean(CPUE,1);
y1e = np.std(CPUE,1);
y2m = np.mean(Tau,1);
y2e = np.std(Tau,1);


## Plot Connectance over evo time
fig,ax1 = plt.subplots(1, figsize=(10, 8),edgecolor=[.4,.4,.4]);
ax1.errorbar(x,y1m,yerr=y1e,xerr=None,\
        alpha=1.,color=[1,.7,.25],lw=2,fmt='-o',ms=10,mec=[1,.7,.25]);
for tl in ax1.get_yticklabels():
	    tl.set_color([1,.7,.25])

ax2 = ax1.twinx()
ax2.axis = [-0.1,1.1,0.012,0.018]
ax2.errorbar(x,y2m,yerr=y2e,xerr=None,\
        alpha=1.,color=[.3,.3,.3],lw=2,fmt='-o',ms=10,mec=[.3,.3,.3]);
for tl in ax2.get_yticklabels():
	    tl.set_color([.3,.3,.3])

ax1.set_xlabel('Cooperation', color='k',fontsize=14)
ax1.set_ylabel('CPUE', color=[1,.7,.25],fontsize=14)
ax2.set_ylabel('Fishing Time', color=[.3,.3,.3],fontsize=14)

plt.savefig('./PDF/Fig_simple.pdf',dpi=600,bbox_inches='tight')

