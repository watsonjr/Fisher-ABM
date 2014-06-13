###### PRAMETERS / CONSTANTS for running the model to equilibrium
function init_equilibrium()

#### Initialize fish
function fnc_init_fish(PF_sig)
    Fish_sx = rand(Float64,PS_n,2).*repmat([GRD_mx GRD_mx],PS_n,1);
    Fish_fx = zeros(Float64,PF_n*PS_n,2); # location of fish
    a = zeros(PS_n,PF_n);
    for i = 1:PS_n
    	a[i,:] = ones(PF_n) .* i;
    end
    Fish_fs = a'[:]; # fish school index
    for i = 1:PF_n*PS_n # place each fish
        Fish_fx[i,:,1] = mod(Fish_sx[Fish_fs[i],:,1] + 
        		randn(1,2).*PF_sig, [GRD_mx GRD_mx]);
    end
    return Fish_fx, Fish_fs, Fish_sx
end
Fish_fx,Fish_fs,Fish_sx = fnc_init_fish(PF_sig);
fish = Fish(Fish_fx,Fish_fs,Fish_sx);

##### Initialize fishers
Cons_H   = zeros(Float64,PC_n); # catch at time t
Cons_xy  = zeros(Float64,PC_n,2); # Fisher locations through time
Cons_xy[:,:,1]=rand(Float64,PC_n,2).*repmat([GRD_mx GRD_mx],PC_n,1);
Cons_s   = randn(PC_n); Cons_s[Cons_s.<0]=-1; Cons_s[Cons_s.>0]=1;
Cons_cn  = zeros(PC_n,PC_n);
Dmin     = zeros(Float64,PC_n); # distance to nearest fish
DXY 	 = [randn(PC_n) randn(PC_n)];
DXY 	 = DXY ./ sqrt(DXY[:,1].^2 + DXY[:,2].^2);
VR  	 = zeros(Float64,PC_n);  # speed modulator (proportional to distance to fish)
RN  	 = zeros(Float64,1,PC_n); # random communication (in fnc_information)
JJ  	 = zeros(Int,PC_n); # index of nearest fish (0 if nothing near)
KK  	 = zeros(Int,PC_n); # index of catch [0,1]
cs  	 = zeros(Float64,PC_n);
D   	 = ones(Float64,PC_n);
cons = Fishers(Cons_xy,Cons_H,Cons_s,Cons_cn,Dmin,DXY,VR,JJ,KK,cs,D);

##### Initialize output storage
OUT        = Output(Fish_fx,Cons_xy,Fish_sx,Cons_H);

return fish,cons,OUT
end

#### REFS
# 1 - http://math.stackexchange.com/questions/52869/numerical-approximation-of-levy-flight/52870#52870
# 2 - http://www.johndcook.com/standard_deviation.html
