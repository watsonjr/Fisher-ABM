

import numpy as np
import igraph
from pylab import *

M = np.load("../Data/Data_WSLS.npz")
SN = M["SOCIAL_NETWORK"]

C  = np.zeros(SN.shape[2])
D = np.zeros(SN.shape[2])
A = np.zeros(SN.shape[2])
M = np.zeros(SN.shape[2])
Fr = np.zeros((SN.shape[0],SN.shape[2]))
for i in np.arange(0,SN.shape[2]):
    g = igraph.Graph.Adjacency(SN[:,:,i].tolist())
    Fr[:,i] = np.sum(SN[:,:,i],1)-1
    C[i] = g.clique_number()
    D[i] = np.mean(g.degree())
    A[i] = g.assortativity_degree()
    #M[i] = g.modularity(g.community_optimal_modularity())
    M[i] = g.modularity(g.community_leading_eigenvector())
    print(i)



##Plot
plt.close('all')
#plt.figure()
#xs, ys = zip(*[(left, count) for left, _, count in 
#    g.degree_distribution().bins()])
#plt.bar(xs, ys)
#plt.title("degree distribution")

Fmax = 10;
FR = np.zeros((Fmax,SN.shape[2]))
for i in np.arange(0,SN.shape[2]):
    fr = np.sum(SN[:,:,i],1)-1
    n,f = np.histogram(fr,np.arange(0,Fmax+1))
    FR[:,i] = n
plt.figure()
plt.pcolormesh(np.log10(FR+1))

figure()
for i in np.arange(Fr.shape[0]):
    plt.plot(Fr[i,:],alpha=0.25)
plt.plot(mean(Fr,0),'r',lw=2)
plt.title("total number of friends over time")

figure()
plt.plot(C)
plt.title("Clique number")

figure()
plt.plot(A)
plt.title("Assortativity")

figure()
plt.plot(D)
plt.title("mean degree")

figure()
plt.plot(mean(Fr,0))
plt.title("mean number of friends")

figure()
plt.plot(M)
plt.title("Optimal modularity")

plt.show()

