
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









