



function ptest(in)
	out = rand() + in
end

@everywhere 
IN  = ones(1000,1000);
OUT = pmap(ptest,IN)


a = randn(1000)
@parallel (+) for i=1:100000
	  f(a[randi(end)])
end


nheads = @parallel (-) for i=1:200000000
	  int(randbool())
end


