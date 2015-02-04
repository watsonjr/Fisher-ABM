#### GENERAL UTILITY FUNCTIONS

macro R_str(s)
    #Small trick for raw strings (e.g. for LaTeX symbols)
    #Write R"\nu" to convert it to "\\nu" which is printed as "\nu"
    s
end

function save_results(res,headers,fname,param=None)
    ## Export results into a text file with name fname, including 
    dim=length(size(res)) #dimension of array
    cellen=length(res[[1 for i in 1:dim]... ] ) #length of cell
    
    if param==None
        param=PRM
    end

    result=zeros(Float64,size(res)... , cellen ) #Printable array containing results
    for i in Iterators.product([1:j for j in size(res)]...) #Multi-index
        result[i...,:]=[x for x in res[i...]]'
    end
    npzwrite("./Data/$(fname).npy", result)
    fout=open("./Data/$(fname).dat", "w")

    write(fout,"# $(headers) \n")
    for f in names(param)
        if ! in("$f",split(headers) )
            write(fout,"$f    $(getfield(PRM,f))\n" )
        end
    end
    close(fout)
end



