#from plot_benchcore import *
import matplotlib.pyplot as plt
import numpy as np

#=========== FIGURE 2A ==================
fig = plt.figure()
ax1 = plt.subplot(1, 1, 1)

#! data
data = np.load("../Data/Data_1fisher.npz")
TS = data["Ts"]
TS_mu = np.mean(TS,1);
xs = data['PC_rp']
TS_an = data['Ts1']

p1 = ax1.plot(xs,TS,'r',alpha=0.15,lw=1,label='realization') 
p2 = ax1.plot(xs,TS_mu,'r',alpha=1,lw=2,label='mean')
p3 = plt.plot(xs,TS_an,'b');
#plt.legend([p1,p2],['1','2'],prop={'size':60}) # legend
#plt.axis([0.1,0.9,0,750])
plt.xlabel(r'$C_p$',fontsize=14,labelpad=10)
plt.ylabel(r'$T_s^*$',fontsize=14,labelpad=-6)
plt.axis([0.1,0.9,0,1000])
plt.show()
fig.savefig('./PNG/Fig_2a.png',dpi=600,bbox_inches='tight')
