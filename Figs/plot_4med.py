#from plot_benchcore import *
#from matplotlib import gridspec
#import mycolor
import numpy as np
import matplotlib.pyplot as plt

#=========== Preamble ==================
#filename="Data_rndcliq_a"
#set_constants(filename)
#data=load_data(filename)
data = np.load("../Data/Data_nfisher_med.npz")

#! stuff to plot
occur=data["occur"] # occurance of a clique size

#clique sizes with more than one occurrence
sizes=np.array([i for i,j in enumerate(occur) if j>1])

#! doing what?
H=data["Hrate_mu"][sizes]
Hvar=data["Hrate_var"][sizes]*occur[sizes]/(occur[sizes]-1) #unbiased variance
Hstd=np.sqrt(Hvar)
TSvar=data["TS_var"][sizes]*occur[sizes]/(occur[sizes]-1) #unbiased variance
TSstd=np.sqrt(TSvar)
TS=data["TS_mu"][sizes]

#! transform TS to encounter rate
tstd_max = TS + TSstd
TS = 1./TS
TSstd = np.abs((1./tstd_max) - TS)

#! plot
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

ax1.errorbar(sizes+1,H,yerr=Hstd,color=[.2,.2,.2],lw=2)
ax2.errorbar(sizes+1,TS,yerr=TSstd,color=[1,.75,.5],lw=2)

i = np.where(H==H.min())
j = np.where(H==H.max())
k = np.where(TS==TS.min())
l = np.where(TS==TS.max())

ax1.axis([1,18,H[i]-Hstd[i],H[j]+Hstd[j]])
ax2.axis([1,18,TS[k]-TSstd[k],TS[l]+TSstd[l]])

ax1.set_xlabel('Clique size',fontsize=13,color=[0,0,0])
ax1.set_ylabel('Average harvest rate',fontsize=13,color=[.2,.2,.2])
ax2.set_ylabel('Average encounter rate',fontsize=13,color=[1,.75,.5])

ax1.spines['left'].set_color([.2,.2,.2])
ax2.spines['right'].set_color([1.,.75,.5])
ax1.tick_params(axis='y',colors=[.2,.2,.2])
ax2.tick_params(axis='y',colors=[1.,.75,.5])

#! save
plt.show()
fig.savefig('./PNG/Fig_4med.png',dpi=600,bbox_inches='tight')



