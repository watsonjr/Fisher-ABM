import numpy as np
import matplotlib.pyplot as plt

## Load in data
Harv = np.load("../Data/Data_harvest.npy")

## Estimate times between hauls
HS = np.empty([Harv.shape[0]],dtype=object)
TT = np.empty([1],dtype=int);
for i in np.arange(Harv.shape[0]):
    a = Harv[i,:]
    b = np.diff(a)
    j = np.where(b==1)[0]
    c = np.diff(j)
    HS[i] = c
    TT = np.concatenate((TT,c))
TT = np.delete(TT,0)

## plot
#plt.hist(np.log10(HS[4]),bins=20)
plt.hist(np.log10(TT),bins=20)
plt.show()


