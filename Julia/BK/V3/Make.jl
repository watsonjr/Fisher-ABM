
#### Preamble
using NPZ, Distance, Devectorize
include("make_parameters.jl");
include("sub_functions.jl");
include("make_season.jl");


#### Run a season
Fish_xy,Cons_xy,Cons_H = make_season(Fish_xy,Cons_xy,Cons_H);


#### Save
npzwrite("./Data/Data_fish.npy", Fish_xy)
npzwrite("./Data/Data_fishers.npy", Cons_xy)
npzwrite("./Data/Data_fclust.npy", Fclust_xy)

#### Useful Julia code
###### Grid fish density
#ii = int(round(Fish_xy[:,:,1]));
#jj,mm,kk = hist2d(ii,Grd_x[1:2:end],Grd_y[1:2:end]);

