###### PARAMETERS / CONSTANTS for running the model to equilibrium

function initialize_e()
    #### Initialize params
    P = T_Params(100.,100.,1.,1.,100.*1.,100.*1.,[0:1.:100.*1.],[0:1.:100.*1.], # Grid
        300,1.,0.01, # Fish
        8,1.5,50.,30.,1.,0.5, # Fishers
        5,[1:5],30.,2.5);  # clusters

    #### Initialize fish
    Fclust_xy = rand(Float64,P.Cl_n,2).*repmat([P.G_mx P.G_my],P.Cl_n,1);
    Fish_xy = zeros(Float64,P.F_n,2); # location of fish
    Fish_cl = zeros(Float64,P.F_n); # index of cluster centre fishers are nearest
    for i = 1:P.F_n # place each fish
        j = ceil(rand(1)*P.Cl_n); # index of which cluster center
        Fish_cl[i] = j[1];
        Fish_xy[i,:,1] = mod(Fclust_xy[j,:,1] + randn(1,2).*P.F_dx, [P.G_mx P.G_my]);
    end
    Fish = T_Fish(Fish_xy,Fish_cl,Fclust_xy);

    ##### Initialize fishers
    Cons_H     = zeros(Int,P.C_n); # catch at time t
    Cons_xy    = zeros(Float64,P.C_n,2); # Fisher locations through time
    Cons_xy[:,:,1]=rand(Float64,P.C_n,2).*repmat([P.G_mx P.G_my],P.C_n,1);
    Con = T_Fishers(Cons_xy,Cons_H);

    #### Define Variables
    Dmin  = zeros(Float64,P.C_n); # distance to nearest fish
    DDx   = zeros(Float64,P.C_n); # x vector compontent to nearest fish
    DDy   = zeros(Float64,P.C_n); # y vector compontent to nearest fish
    ANG   = zeros(Float64,P.C_n); # fisher direction angle
    VR    = zeros(Float64,P.C_n);  # speed modulator (proportional to distance to fish)
    RN    = zeros(Float64,1,P.C_n); # random communication (in fnc_information)
    JJ    = zeros(Int,P.C_n); # index of nearest fish (0 if nothing near)
    KK    = zeros(Int,P.C_n); # index of catch [0,1]
    Var  = T_Vars(Dmin,DDx,DDy,ANG,VR,RN,JJ,KK);

    ##### Initialize waiting times
    Tau_n    = zeros(Int,P.C_n); # number of events (hauls) so far
    Tau_t    = zeros(Int,P.C_n); # time between current events (waiting time)
    Tau_s    = zeros(Float64,P.C_n); # cumulative waiting time
    Tau_mu   = zeros(Float64,P.C_n); # running average waiting time (wait_sum / wait_n)
    Tau_dmu=ones(Float64,P.C_n)*999; # min diff in mean waiting times between t and t-1
    Tau = T_Tau(Tau_n,Tau_t,Tau_s,Tau_mu,Tau_dmu);

    return P,Fish,Con,Var,Tau
end

#### REFS
# 1 - http://math.stackexchange.com/questions/52869/numerical-approximation-of-levy-flight/52870#52870
