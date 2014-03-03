import numpy as np
import matplotlib.pyplot as plt

## Load in data
Harv = np.load("../Data/Data_harvest.npy")

## Data to plot
y = Harv[2,:];
x = np.arange(0,y.size);
d = y*x;
d = d[d>0];
bins = np.arange(0,y.size+1,20); ##<< bin width
hist_1, edges = np.histogram(d,bins);
cen_1 = bins[1:] - (np.diff(bins)/2);


## Calculate times between catches
d = np.hstack((0,d))
tau = np.diff(d);
bins = np.linspace(0,tau.max(),num=20);
hist_2, edges = np.histogram(d,bins); # time between catches
cen_2 = bins[1:] - (np.diff(bins)/2);


## plot harvest stems
fig = plt.figure(1, figsize=(10, 8),edgecolor=[.4,.4,.4]);
markerline, stemlines, baseline = plt.stem(cen_1[hist_1>0],hist_1[hist_1>0], '-');
plt.setp(stemlines,color='g',lw=1.5);
plt.setp(markerline,mfc='b',mec='b',markersize=10);
plt.savefig('./PNG/Fig_harvest.png',dpi=100,bbox_inches='tight')
#plt.close('all')

## plot time between catch
fig = plt.figure(1, figsize=(10, 8),edgecolor=[.4,.4,.4]);
markerline, stemlines, baseline = plt.stem(cen_2,hist_2, '-');
plt.setp(stemlines,color='g',lw=1.5);
plt.setp(markerline,mfc='b',mec='b',markersize=10);
plt.savefig('./PNG/Fig_tau.png',dpi=100,bbox_inches='tight')
#plt.close('all')


