###### PRAMETERS / CONSTANTS for running the model to equilibrium
function init_equilibrium()

#### Initialize fish
#! all we need are fish location (x,y)
#! the index of the school each fish is attached to (ind)
#! and the school location (x,y)
#! and the KDTREE
fish_sx = rand(Float64,PS_n,2).*repmat([GRD_mx GRD_mx],PS_n,1);
fish_fx = zeros(Float64,PF_n*PS_n,2); # location of fish
a = zeros(PS_n,PF_n);
for i = 1:PS_n; a[i,:] = ones(PF_n) .* i; end
fish_fs = a'[:]; # fish school index
for i = 1:PF_n*PS_n; fish_fx[i,:,1] = mod(fish_sx[fish_fs[i],:,1] + 
			randn(1,2).*PF_sig, [GRD_mx GRD_mx]); end

##### Initialize fishers
#! all we need are the locations of each fisher (x,y)
#! the index of the nearest fish (from KDTree)
#! the distance to that nearest fish
#! the unit vector pointing at that fish (random at start)
#! the harvest count (1/0)
#! the social strategy (make/break friendships)
#! an index of whether, in the absense of fish, they steam or search
#! the social network
cons_xy  = rand(Float64,PC_n,2).*repmat([GRD_mx GRD_mx],PC_n,1);
cons_Ni  = Array(Int,PC_n,1);
cons_Nd  = Array(Float64,PC_n,1)
DXY      = [randn(PC_n) randn(PC_n)];
cons_dx  = DXY ./ sqrt(DXY[:,1].^2 + DXY[:,2].^2);
cons_H   = zeros(Float64,PC_n); # catch at time t
cons_s   = randn(PC_n);cons_s[cons_s.<0]=-1;cons_s[cons_s.>0]=1;cons_s=int(cons_s)
cons_mi  = int(zeros(PC_n));
cons_sn  = ones(PC_n,PC_n) .* eps(); 
for j = 1:PC_n; cons_sn[j,j] = 1; end;


##### Initialize
fish = Fish(fish_fx,fish_fs,fish_sx);
cons = Fishers(cons_xy,cons_Ni,cons_Nd,cons_dx,cons_H,cons_s,cons_mi,cons_sn)
OUT  = Output(fish_fx,cons_xy,fish_sx,cons_H);

return fish,cons,OUT
end

#### REFS
# 1 - http://math.stackexchange.com/questions/52869/numerical-approximation-of-levy-flight/52870#52870
# 2 - http://www.johndcook.com/standard_deviation.html
