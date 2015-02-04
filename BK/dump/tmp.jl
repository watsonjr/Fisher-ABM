

using NearestNeighbors
using Distance_periodic

## Naive
X = [1 1;5 5;98 98]';
v = X[:, 1]
t1 = NaiveNeighborTree(X,Euclidean())
inds, dists = nearest(t1, v, 1)



