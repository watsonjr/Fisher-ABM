from plot_benchcore import *
from matplotlib import gridspec
import matplotlib.colors as mcolors
import mycolor

#=========== CUSTOM COLORMAP =========#
levs = range(24)
assert len(levs) % 2 == 0, 'N levels must be even.'
cmap1 = mcolors.LinearSegmentedColormap.from_list(name='red_white_blue', 
                                                 colors =[(0, 0, 1), 
                                                          (1, 1., 1), 
                                                          (1, 0, 0)],
                                                 N=len(levs)-1,
                                                 )
cmap2 = mcolors.LinearSegmentedColormap.from_list(name='white_blue', 
                                                 colors =[(0., 0., 1.), 
                                                          (1, 1, 1)],
                                                 N=len(levs)-1,
                                                 )

#=========== Preamble ==================
#fig = plt.figure(figsize=(10,5))
filename="Data_2fisher_1school_inf"
set_constants(filename)
data_inf=load_data(filename)
filename="Data_2fisher_1school_noinf"
data_noinf=load_data(filename)

TS = data_noinf['\\tau_s^R']
TSinf = data_inf['\\tau_s^R']
H = data_noinf['H']
Hinf = data_inf['H']
X0 = data_inf['PS_p']
Y0 = data_inf['PC_q']

X = 1./X0
Y = constants["PF_n"]/Y0
taus1 = X.copy() # optimal for 1 fisher
for i in range(taus1.shape[0]):
    for j in range(taus1.shape[1]):
        args=get_constants(PS_p=data_inf['PS_p'][i,j],PC_q=data_inf['PC_q'][i,j],   PC_rp=data_inf['PC_rp'][i,j] )
        taus1[i,j]=tausr(*args)
X1,Y1 = np.log10(X/taus1),np.log10(Y/taus1)

###### Find maxima and minima in VOI for mfisher_nclique exp
Z = Hinf/H
i,j = np.where(Z == Z.max());
PS_p_max = X0[i,j] #PS_p
PC_q_max = Y0[i,j] #PC_q
X1[i,j][0]
Y1[i,j][0]
i,j = np.where(Z == Z.min());
PS_p_min = X0[i,j] #PS_p
PC_q_min = Y0[i,j] #PC_q

###### Heatmap for H
fig = plt.figure()
Z = Hinf/H
plt.xlabel(r"$log_{10}(\tau_l/\tau_s^1)$",fontsize=15)
plt.ylabel(r"$log_{10}(\tau_h/\tau_s^1)$",fontsize=15)
cmap = plt.get_cmap('bwr')
plt.pcolormesh(X1, Y1, Z,  vmin=(Z).min(), vmax=(Z).max(),cmap="seismic")
plt.clim([0.8,1.2]) # make it symmetric
plt.colorbar()
plt.xlim(X1.min(),X1.max()); plt.ylim(Y1.min(),Y1.max())
plt.show()
fig.savefig('./EPS/Fig_3a.eps',dpi=600,bbox_inches='tight')

###### Heatmap for Ts
fig = plt.figure()
Z = TSinf/TS
plt.xlabel(r"$log_{10}(\tau_l/\tau_s^1)$",fontsize=15)
plt.ylabel(r"$log_{10}(\tau_h/\tau_s^1)$",fontsize=15)
cmap = plt.get_cmap('Blues_r')
plt.pcolormesh(X1, Y1, Z,  vmin=(Z).min(), vmax=(Z).max(), cmap="seismic")
plt.clim([.8,1.2])
plt.colorbar()
plt.xlim(X1.min(),X1.max()); plt.ylim(Y1.min(),Y1.max())
plt.show()
fig.savefig('./EPS/Fig_3b.eps',dpi=600,bbox_inches='tight')



