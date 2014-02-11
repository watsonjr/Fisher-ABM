###### Parameters for fish agent-based model

#### Domain parameters
Grd_nx = 100.;
Grd_ny = 100.;
Grd_dx = 1.;
Grd_dy = 1.;
Grd_mx = Grd_nx.*Grd_dx;
Grd_my = Grd_ny.*Grd_dy;
Grd_x  = [0:Grd_dx:Grd_mx];
Grd_y  = [0:Grd_dy:Grd_my];

#### Model parameters
P_Tend = 1000; # number of timesteps
P_ddx  = 10.; # P_ddx/2 = the dx and dy of fish finder
P_catch = 4; # distance at which fishers can catch fish
P_prob  = 0.5; # probability of catching fish once near

#### Fish parameters
## cluster movement - see ref 1
# theta=rand*2*pi; # direction to move
# f=rand^(-1/alpha); # step size
# x(n)=x(n-1)+f*cos(theta); # update x
# y(n)=y(n-1)+f*sin(theta); # update y
Fc_n   = 5;    # number of cluster centres
Fc_ang = 30;   # cluster angle change per unit time (i.e. +- Fc_ang)
Fc_alp = 1.5;  # levy flight power law constant (determins max dist moved per unit time)
Fc_min = .01;    # min distance that clusters can move per unit time

Fclust = zeros(Fc_n,2,P_Tend);
Fclust[:,:,1] = rand(Float64,Fc_n,2).*repmat([Grd_nx*Grd_dx Grd_ny*Grd_dy],Fc_n,1);

## fish
F_n    = 500;  # number of fish
F_dx   = 2;  # distance from cluster center
F_id   = zeros(F_n,1);

## Initial fish locations
Fish = zeros(F_n,2,P_Tend);
for i = 1:F_n # place each fish
    j = ceil(rand(1)*Fc_n); # index of which cluster center
    Fish[i,:,1] = mod(Fclust[j,:,1] + randn(1,2).*F_dx, [Grd_mx Grd_my]);
end

#### Fisher parameters
C_n     = 5; # number of fishers
C_v     = .5; # max speed
C_dx    = 5; # width of fish finder

## Initialize fishers
Cons    = zeros(C_n,2,P_Tend);
Cons[:,:,1] = rand(Float64,C_n,2).*repmat([Grd_nx*Grd_dx Grd_ny*Grd_dy],C_n,1);
Harv    = zeros(C_n,P_Tend);

## Stuff I use
DXY = zeros(C_n,2);
Cnr = zeros(C_n);
theta_f = zeros(F_n,1);
theta_c = zeros(C_n,1);

#### REFS
# 1 - http://math.stackexchange.com/questions/52869/numerical-approximation-of-levy-flight/52870#52870
