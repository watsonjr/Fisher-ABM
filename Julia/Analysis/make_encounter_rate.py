import numpy as np
import matplotlib.pyplot as plt

## Load in data
Harv = np.load("../Data/Data_fisher_H.npy")

## Calculate times between catches
MU = np.zeros(Harv.shape[0]);
MU_s = np.zeros(Harv.shape[0]);
for j in np.arange(0,Harv.shape[0]):

    y = Harv[j,:];
    x = np.arange(0,y.size);
    d = y*x;
    d = d[d>0];
    d = np.hstack((0,d))
    tau = np.diff(d);

    ## Calculate sampling mean and error
    mu = np.zeros(tau.size); sd = np.zeros(tau.size);
    for i in np.arange(0,tau.size):
        mu[i] = np.mean(tau[0:i]);
        sd[i] = np.std(tau[0:i]);

    ## Perform running sampling mean
    su = np.zeros(tau.size); mu_s = np.zeros(tau.size);
    for i in np.arange(1,tau.size):
        su[i]   = su[i-1] + tau[i-1];
        mu_s[i] = su[i] / (i+1);

    MU[j] = mu[-1];
    MU_s[j] = mu_s[-1];



## Plot
#  remember: stationary = constant mean, constant variance
fig = plt.figure(1, figsize=(10, 8),edgecolor=[.4,.4,.4]);
plt.plot((mu[10:]));
plt.plot((mu_s[10:]));

fig = plt.figure(2, figsize=(10, 8),edgecolor=[.4,.4,.4]);
plt.plot(np.diff(sd[10:]));

fig = plt.figure(2, figsize=(10, 8),edgecolor=[.4,.4,.4]);
markerline, stemlines, baseline = plt.stem(cen_2,hist_2, '-');
plt.setp(stemlines,color='g',lw=1.5);
plt.setp(markerline,mfc='b',mec='b',markersize=10);
