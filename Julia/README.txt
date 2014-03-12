

#### Useful Julia code
###### Grid fish density
#ii = int(round(Fish_xy[:,:,1]));
#jj,mm,kk = hist2d(ii,Grd_x[1:2:end],Grd_y[1:2:end]);

#(Tau_n,Tau_t,Tau_s,Tau_mu,Tau_dmu,
#   Dmin,DDx,DDy,ANG,VR,RN,JJ,KK,
#   Fish_xy,Fish_ci,Fish_cl,
#   Cons_xy,Cons_H) = init_equilibrium();
#
#fish = Fish(copy(Fish_xy),copy(Fish_ci),copy(Fish_cl));
#cons = Fishers(copy(Cons_xy),copy(Cons_H));
#vars=Vars(copy(Dmin),copy(DDx),copy(DDy),copy(ANG),copy(VR),copy(RN),copy(JJ),copy(KK));
#tau = Tau(copy(Tau_n),copy(Tau_t),copy(Tau_s),copy(Tau_mu),copy(Tau_dmu));

#npzwrite("./Data/Data_fclust.npy", fish.cl)
#npzwrite("./Data/Data_fisher_H.npy", cons.H)


