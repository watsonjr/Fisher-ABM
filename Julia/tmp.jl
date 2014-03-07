
using PyPlot, NPZ
SN = npzread("./Data/Data_SN.npy");

C = sum(SN,1);
C = sum(C,2);

D = sum(SN,2);
D = mean(D,1);

c = zeros(size(SN,3));
d = zeros(size(SN,3));
c[:,1] = C[1,1,:];
d[:,1] = D[1,1,:];

plot(c);
plot(d);



