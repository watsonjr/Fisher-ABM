

##########################!!!!! SIMULATION EXPERIMENTS !!!!!!!!################


##################### BASIC TESTS
    
function do_timingtest()
    global PC_rp = 0.92; # choose random change in walk
    school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
    @time  make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
    npzwrite("./Data/Data_fish.npy", OUT.fish_xy)
    npzwrite("./Data/Data_fishers.npy", OUT.cons_xy)
    npzwrite("./Data/Data_clusters.npy", OUT.schl_xy)
    npzwrite("./Data/Data_harvest.npy", OUT.cons_H)
end


##################### 1 SCHOOl 1 FISHER
    
function do_first_passage()
    ###### FIRST PASSAGE time (stop as soon as a school is found)
    println("First passage time")
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
        Ts[i] = mean(time),std(time);#cons.measure["Ts"], cons.measure["Tv"];
        print(i,"\n")
    end

    TS = zeros(Float64,size(RP,1),PC_n)
    for i = 1:length(RP)
        TS[i,:] = Ts[i][1]
    end
    npzwrite("./Data/Data_firstpass.npy", TS)
    npzwrite("./Data/Data_firstpass_xs.npy", RP)
end

function do_fig2a()
    ######  TYPICAL SEARCH time over many searches
    println("Typical search time")
    RP = linspace(0.2,.9,20);
    result = cell(size(RP));

    for i = 1:length(RP)
        global PC_rp = RP[i];
        school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
        #FLAGS["spying"]=true
        make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
        result[i] = cons.measure["Ts"], cons.measure["Tv"],cons.measure["f1"];
        print(i,"\n")
    end

    TS = zeros(Float64,size(RP,1),PC_n)
    F1 = zeros(Float64,size(RP,1),PC_n)
    for i = 1:length(RP)
        TS[i,:] = result[i][1]
        F1[i,:] = result[i][3]
    end
    npzwrite("./Data/Data_Fig2a.npy", TS)
    npzwrite("./Data/Data_Fig2a_f1.npy", F1)
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
            Ts[i,j] =[SIG[i],FF[j]], cons.measure["Ts"], cons.measure["Tv"],PC_rp;
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


##################### 1 SCHOOL 2 FISHERS
######## Test effect of friendship, on Tau_s and CPUE, as Tau_l is varied

function do_fig3()
    ###### Test performance as a function of tau_l and tau_h
    ## basic values of C_f and F_sig
    global PC_f=GRD_nx*0.05
    global PF_sig=GRD_nx*0.05
    
    JUMP = logspace(log10(0.001),log10(.1),6);
    CATCH  = logspace(log10(0.01),log10(1),6);
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
            result1=cons.measure["Ts"]

            println("With full information")
            #With information
            global PC_lambda=1;
            school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
#            println(cons.SN)
            make_season(school,fish,cons,fishtree,EVENTS,FLAGS);
            result2=cons.measure["Ts"]

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


######## Test effect of friendship, on Tau_s and CPUE, as F_n is varied

##################### 1 FISH N FISHERS (do 5, 10, 20 fishers)
######## Stupid Optimal Social Network test (just keep adding links)

######## Evolving social network - find the ESS 








######## END  #######





## simulate a simple scenario
function sim_simple()
sn = linspace(1e-6,1,10);   # types of prosociality
trips = 20; # number of repeats
CPUE = Array(Float64,length(sn),trips);
Tau  = Array(Float64,length(sn),trips);
for i = 1:length(sn)
    ## modulate social network
    SN = ones(PC_n,PC_n) .* sn[i];
    for k = 1:PC_n; SN[k,k] = 1; end;

    for j = 1:trips

        ## run model
        school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
        cons.SN=SN
        FLAGS["save"]=true
        make_season(school,fish,cons,fishtree,EVENTS,FLAGS,1);

        ## record
        CPUE[i,j] = mean(cons.cs ./ cons.Dist);
        Tau[i,j]  = mean(cons.Dist);
    end
    print(i/length(sn))

end
return CPUE, Tau
end



## simulate a greedy search for the social network that maximizes the 
#! FLEET's average catch per unit effort
function sim_fleet()

seasons = 1
trips = 12; # number of repeats
Social_network = Array(Float64,PC_n,PC_n,seasons)
cpue = Array(Float64,trips);
CPUE = Array(Float64,seasons); CPUE[1] = 0;
STR = 1; # initial strategy (1=makefriends,0=breakfriends)

for i = 2:seasons # greedy search over seasons
    
    ## Update social network according to strategy
    #! maybe break symmetry??
    if STR == 1 # make friends
        k = find(SN.==eps());
        k = k[ceil(rand().*length(k))];
        SN[k] = 1;
        #k = ind2sub(size(SN),k);
        #SN[k[1],k[2]] = 1;
        #SN[k[2],k[1]] = 1;
    elseif STR == 0 #break friends
        k = find(SN.==1);
        k = k[ceil(rand().*length(k))];
        k = ind2sub(size(SN),k);
        if k[1]!=k[2]
            SN[k[1],k[2]] = eps();
            #SN[k[2],k[1]] = eps();
        end
    end
    Social_network[:,:,i] = SN;

    ## Run fishing seasons
    for j = 1:trips # build up catch statistics

        ## run model
        school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
        cons.SN=SN
        make_season(school,fish,cons,fishtree,EVENTS,FLAGS);

        ## record
        cpue[j] = mean(cons.cs ./ cons.Dist);
    end

    ## Update fleet performance and strategy - win stay, loose shift
    CPUE[i] = mean(cpue);
    if (CPUE[i]-CPUE[i-1]) > 0 # if CPUE increased
        STR = STR; # stay 
    else # else if CPUE decreased
        STR = abs(STR-1); # change
    end

    ## Ticker
    #print(i/seasons,"\n")
end
return CPUE,Social_network
end


## simulate a greedy search for the social network that maximizes the 
#! INDIVIDUAL's average catch per unit effort
function sim_individual()

seasons = 200
trips   = 30; # number of repeats
Social_network = Array(Float64,PC_n,PC_n,seasons)
cpue = Array(Float64,PC_n,trips);
CPUE = Array(Float64,PC_n,seasons); CPUE[:,1] = 0;
STR  = ones(PC_n); # initial strategy (1=makefriends,0=breakfriends)

for T = 2:seasons # greedy search over seasons

    ## Update social network according to strategy
    #! find all fishers who want to make friends
    #! pair them up
    #! those who want to break friendships, just do it randomly

    #! MAKE friends
    idy   = find(STR.==1) # those fishers who want to make friends
    if isempty(idy) == 0
        LL       = length(idy);
        idx   = randperm(LL); # randomly associate pairs to become friends
        if mod(LL,2) == 0
            idx = reshape(idx,(int(LL/2),2)); # if even number
        else
            idx=idx[1:end-1]; # if odd number 
            idx = reshape(idx,(int((LL-1)/2),2));
        end
        for i = 1:size(idx,1);
            SN[idx[i,1],idx[i,2]] = 1.;
            SN[idx[i,2],idx[i,1]] = 1.;
        end
    end

    #! BREAK friends
    idy   = find(STR.==0) # those fishers who want to break friendships
    if isempty(idy) == 0
        LL       = length(idy);
        for i = 1:LL
            j = find(SN[idy[i],:] .== 1.); # find your current friends
            j = j[j.!=i]; # don't break up with yourself
            if isempty(j)==0
                j = j[randperm(length(j))[1]] # choose one at random
                SN[i,j] = eps();
                SN[j,i] = eps();
            end
        end
    end


    ## Run fishing seasons
    for j = 1:trips # build up catch statistics

        ## run model
        school,fish,cons,fishtree,EVENTS,FLAGS,OUT = init_equilibrium();
        cons.SN=SN
        make_season(school,fish,cons,fishtree,EVENTS,FLAGS);

        ## record
        cpue[:,j] = cons.cs ./ cons.Dist;
    end

    ## Update fleet performance and strategy - win stay, loose shift
    Social_network[:,:,T] = SN;
    CPUE[:,T] = mean(cpue,2);
    for i = 1:PC_n
        if (CPUE[i,T]-CPUE[i,T-1]) > 0 # if CPUE increased
            STR[i] = STR[i]; # stay
        else # else if CPUE decreased
            STR[i] = abs(STR[i]-1); # change
        end
    end

    ## Ticker
    print(T/seasons,"\n")
end

return CPUE,Social_network
end









