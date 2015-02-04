#from matplotlib import gridspec
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

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
data = np.load("../Data/Data_2fisher.npz")
TS_noinf = 1./data['TS_noinf']
TS_inf = 1./data['TS_inf']
H_noinf = data['H_noinf']
H_inf = data['H_inf']
taul = data['taul']
tauh = data['tauh']
taus1 = data['taus1']

X = np.log10(taul / taus1)
Y = np.log10(tauh / taus1)

###### Find location of quantiles
Z = H_inf/H_noinf
i_max,j_max = np.where(Z==Z.max())
i_min,j_min = np.where(Z==Z.min())
dz = Z - np.percentile(Z.flatten(),60)
i_med,j_med = np.where(np.abs(dz) == np.abs(dz).min())

dx = np.abs(np.diff(X[:,0])[0])/2
dy = np.abs(np.diff(Y[0,:])[0])/2

###### Heatmap for H
fig = plt.figure()
Z = H_inf/H_noinf
plt.xlabel(r"$log_{10}(\tau_l/\tau_s^1)$",fontsize=15)
plt.ylabel(r"$log_{10}(\tau_h/\tau_s^1)$",fontsize=15)
cmap = plt.get_cmap('bwr')
plt.pcolormesh(X,Y, Z,vmin=(Z).min(), vmax=(Z).max(),cmap="seismic")
plt.clim([0.7,1.3]) # make it symmetric
plt.colorbar()

plt.scatter(X[i_max[0],j_max[0]]-dx,Y[i_max[0],j_max[0]]-dy,s=100,
            facecolors='none',edgecolors='w',linewidth=3,alpha=1.)
plt.scatter(X[i_min[0],j_min[0]]-dx,Y[i_min[0],j_min[0]]-dy,s=100,
            facecolors='none',edgecolors=[0,1,0],linewidth=3,alpha=1.)
plt.scatter(X[i_med[0],j_med[0]]-dx,Y[i_med[0],j_med[0]]-dy,s=100,
            facecolors='none',edgecolors='m',linewidth=3,alpha=1.)

plt.xlim(X.min(),X.max()); plt.ylim(Y.min(),Y.max())
plt.show()
fig.savefig('./EPS/Fig_3a.eps',dpi=600,bbox_inches='tight')
fig.savefig('./PNG/Fig_3a.png',dpi=600,bbox_inches='tight')

###### Heatmap for Ts
fig = plt.figure()
Z = TS_inf/TS_noinf
plt.xlabel(r"$log_{10}(\tau_l/\tau_s^1)$",fontsize=15)
plt.ylabel(r"$log_{10}(\tau_h/\tau_s^1)$",fontsize=15)
cmap = plt.get_cmap('Blues_r')
plt.pcolormesh(X, Y, Z,  vmin=(Z).min(), vmax=(Z).max(), cmap="seismic")
plt.clim([.7,1.3])
plt.colorbar()
plt.xlim(X.min(),X.max()); plt.ylim(Y.min(),Y.max())
plt.show()
fig.savefig('./EPS/Fig_3b.eps',dpi=600,bbox_inches='tight')
fig.savefig('./PNG/Fig_3b.png',dpi=600,bbox_inches='tight')



