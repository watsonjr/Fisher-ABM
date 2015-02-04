addprocs(2)
require("tmp.jl")

#! serial
N = 100;
Y = zeros(N)
for i = 1:N
	Y[i] = fnc(i)
end

#! parfor
@parallel for i = 1:N
	j = i
	Y[j] = fnc(j)
end

#! pmap
Y = pmap(fnc,[1:N])



B=1000
samples_pmap = pmap(bootstrapSamples,[b,b,b,b])
samples = vcat(samples_pmap[1],samples_pmap[2],samples_pmap[3],samples_pmap[4])


