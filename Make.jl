    
#### Add modules
using PyPlot, NPZ, Devectorize
using Types, Constants
using PyCall
@pyimport rtree.index as pyrtree

#### Add functions and routines
include("sub_init.jl");
include("sub_functions.jl");
include("sub_routines.jl");
include("Experiments.jl");

##########################!!!!! SIMULATION EXPERIMENTS !!!!!!!!################

#Switches for various experiments

timingtest=false
fig2a=true
fig2b=false
fig3=false

###### run for one season

function do_timingtest()
    global PC_rp = 0.9; # choose random change in walk
    school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
    @time  make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
    npzwrite("./Data/Data_fish.npy", OUT.fish_xy)
    npzwrite("./Data/Data_fishers.npy", OUT.cons_xy)
    npzwrite("./Data/Data_clusters.npy", OUT.schl_xy)
    npzwrite("./Data/Data_harvest.npy", OUT.cons_H)
end


function alt_do_fig2a()
    ##################### 1 FISH 1 FISHER
    ###### Test performance as a function of random walk
    RP = linspace(0.1,.95,20);
    Ts = cell(size(RP));
    global PS_p=0.0
    for i = 1:length(RP)
        global PC_rp = RP[i];
        time=[]
        for iter in 1:300
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            time=vcat(time, make_season(school,fish,cons,fishtree,EVENTS,FLAGS,3) ) ;
        end
        Ts[i] = mean(time),std(time);#cons.Ts, cons.Tv;
        print(i,"\n")
    end

    TS = zeros(Float64,size(RP,1),PC_n)
    for i = 1:length(RP)
        TS[i,:] = Ts[i][1]
    end
    npzwrite("./Data/Data_Fig2a.npy", TS)
    npzwrite("./Data/Data_Fig2a_xs.npy", RP)
end


function do_fig2a()
    ##################### 1 FISH 1 FISHER
    ###### Test performance as a function of random walk
    RP = linspace(0.2,.9,20);
    Ts = cell(size(RP));

    for i = 1:length(RP)
        global PC_rp = RP[i];
        school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
        #FLAGS["spying"]=true
        make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
        Ts[i] = cons.Ts, cons.Tv;
        print(i,"\n")
    end

    TS = zeros(Float64,size(RP,1),PC_n)
    for i = 1:length(RP)
        TS[i,:] = Ts[i][1]
    end
    npzwrite("./Data/Data_Fig2a.npy", TS)
    npzwrite("./Data/Data_Fig2a_xs.npy", RP)
end

function do_fig2b()
    ###### Test performance as a function of C_f and F_sig
    SIG = linspace(1,GRD_mx/10,30);
    FF  = linspace(1,GRD_mx/10,30);
    Ts = cell(size(SIG,1),size(FF,1));
    for i = 1:length(SIG)
        for j = 1:length(FF)
            global PF_sig = SIG[i];
            global PC_f   = FF[j];
            global PC_rp = fnc_optimal_PCrp();
            print(i," ",j,"\n")
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            Ts[i,j] =[SIG[i],FF[j]], cons.Ts, cons.Tv,PC_rp;
        end
    end
    TS = zeros(Float64,size(SIG,1),size(FF,1),PC_n)
    xs = zeros(Float64,size(SIG,1),size(FF,1),2)
    for i = 1:size(TS,1)
        for j=1:size(TS,2)
            TS[i,j,:] = Ts[i,j][2]
            xs[i,j,:]=Ts[i,j][1]
        end
    end
    npzwrite("./Data/Data_Fig2b.npy", TS)
    npzwrite("./Data/Data_Fig2b_xs.npy", xs)
end


function do_fig3()
    ###### Test performance as a function of tau_l and tau_h
    ## basic values of C_f and F_sig
    global PC_f=GRD_nx*0.05
    global PF_sig=GRD_nx*0.05
    
    JUMP = logspace(log10(0.001),log10(.1),30);
    CATCH  = logspace(log10(0.01),log10(1),30);
    Ts = cell(size(JUMP,1),size(CATCH,1));
    for i = 1:length(JUMP)
        for j = 1:length(CATCH)
            global PS_p = JUMP[i];
            global PC_q = CATCH[j];
            global PC_rp = fnc_optimal_PCrp();
            print("$i $j, $PS_p $PC_q \n")
            
            println("Without information")
            #Without information
            global PC_lambda=0;
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
#            println(cons.SN)
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            result1=cons.Ts

            println("With information")
            #With information
            global PC_lambda=1;
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
#            println(cons.SN)
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            result2=cons.Ts

            Ts[i,j,:] = [PS_p,PC_q],result1,result2;
        end
    end
    TS = zeros(Float64,size(JUMP,1),size(CATCH,1),PC_n)
    TSinf = zeros(Float64,size(JUMP,1),size(CATCH,1),PC_n)
    xs = zeros(Float64,size(JUMP,1),size(CATCH,1),2)
    for i = 1:size(TS,1)
        for j=1:size(TS,2)
            TS[i,j,:] = Ts[i,j][2]
            TSinf[i,j,:] = Ts[i,j][3]
            xs[i,j,:]= Ts[i,j][1]
        end
    end
    npzwrite("./Data/Data_Fig3_noinf.npy", TS)
    npzwrite("./Data/Data_Fig3_inf.npy", TSinf)
    npzwrite("./Data/Data_Fig3_xs.npy", xs)
end

if timingtest
    do_timingtest()
end
if fig2a
    do_fig2a()
end
if fig2b
    do_fig2b()
end
if fig3
    do_fig3()
end

##################### 1 FISH 2 FISHERS
######## Test effect of friendship, on Tau_s and CPUE, as Tau_l is varied


######## Test effect of friendship, on Tau_s and CPUE, as F_n is varied






##################### 1 FISH N FISHERS (do 5, 10, 20 fishers)
######## Stupid Optimal Social Network test (just keep adding links)

######## Evolving social network - find the ESS 








####### Simple problem (raise SN entirely)
#@time CPUE,Tau = sim_simple();
#npzwrite("./Data/Data_simple.npz", ["x"=>1, "CPUE"=>CPUE, "Tau"=>Tau])
#
####### Optimize SN for the fleet catch
#@time (CPUE,Social_network) = sim_fleet();
#
####### Optimize an Individual's social network
#@time (CPUE,Social_network) = sim_individual()
#npzwrite("./Data/Data_individual.npz", ["x"=>1, "CPUE"=>CPUE, "SN"=>Social_network])
######## END  #######


