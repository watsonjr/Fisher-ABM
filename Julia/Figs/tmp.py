
OUT = zeros(200)
for i in np.arange(0,199):
    OUT[i] = nx.average_clustering(nx.Graph(SN[:,:,i]))


