from plot_benchcore import *
from matplotlib import gridspec
import matplotlib 
import mycolor

#=========== FIGURE 2A ==================
fig = plt.figure()
ax1 = plt.subplot(1, 1, 1)
filename="Data_1fisher_1school_a"
set_constants(filename)
data=load_data(filename) 
TS=data['\\tau_s^R']
TS_mu = np.mean(TS,1);
xs=data['PC_rp']
TS_an = np.asarray([tausr(*get_constants(PC_rp=x[0])) for x in xs])

p1 = ax1.plot(xs,TS,'r',alpha=0.15,lw=1,label='realization') 
p2 = ax1.plot(xs,TS_mu,'r',alpha=1,lw=2,label='mean')
p3 = plt.plot(xs,TS_an,'b');
#plt.legend([p1,p2],['1','2'],prop={'size':60}) # legend
#plt.axis([0.1,0.9,0,750])
plt.xlabel(r'$C_p$',fontsize=14,labelpad=10)
plt.ylabel(r'$T_s^*$',fontsize=14,labelpad=-6)
plt.show()
fig.savefig('./EPS/Fig_2a.eps',dpi=600,bbox_inches='tight')


##=========== FIGURE 2B ==================
#fig = plt.figure()
#ax1 = plt.subplot(1, 1, 1)
#filename="Data_1fisher_1school_b"
#set_constants(filename)
#data=load_data(filename) 
#GRD_mx = constants["GRD_mx"]
#
#X=data['PF_sig']
#Y=data['PC_f']
#
#TS=data['\\tau_s^R']
#Cp=data['PC_rp']
#Z = np.log10(TS)
##Z = Cp;
#
#plt.pcolormesh(X/GRD_mx,Y/GRD_mx,Z,vmin=abs(Z).min(),vmax=abs(Z).max())
#plt.xlim([0.02,0.1])
#plt.ylim([0.02,0.1])
#plt.xlabel(r'$F_{\sigma} / X$',fontsize=14,labelpad=10)
#plt.ylabel(r'$C_f / X$',fontsize=14,labelpad=-6)
#
#plt.colorbar(label=r'T_s^1')
#ax1.set_rasterized('True')
#plt.show()
#fig.savefig('./EPS/Fig_2b.eps',dpi=600,bbox_inches='tight')
#
