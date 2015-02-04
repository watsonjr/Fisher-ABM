

abstract PreMetric
abstract SemiMetric <: PreMetric
abstract Metric <: SemiMetric
result_type(::PreMetric, T1::Type, T2::Type) = Float64

type Euclidean <: Metric end

#function euclidean(a::AbstractVector,b::AbstractVector)
function evaluate(dist::Euclidean, a::AbstractVector, b::AbstractVector)
	n = length(a)::Int
	s = 0.

	for i = 1:n
		
		## distances
		Dx = a[i] - b[i]

		## periodic boundaries
		c = 50. # constant = half grid dim
		if (abs(Dx).>c) == 1
			Dx = -sign(Dx)*((c*2) - abs(Dx))
		end

		## sum sqr distance
		@inbounds s += (abs(Dx)).^2
	end
	D = sqrt(s)
end
#cityblock(a::AbstractVector, b::AbstractVector) = evaluate(Cityblock(), a, b)
#evaluate(dist::Euclidean,a::AbstractVector,b::AbstractVector)= euclidean(a,b)
euclidean(a::AbstractVector,b::AbstractVector) = evaluate(Euclidean(),a,b)


#function sumsqdiff(a::AbstractVector, b::AbstractVector)
#    n = length(a)::Int
#    s = 0.
#    for i = 1:n
#        @inbounds s += abs2(a[i] - b[i])
#    end
#    s
#end
#evaluate(dist::Euclidean, a::AbstractVector, b::AbstractVector) = sqrt(sumsqdiff(a, b))
#euclidean(a::AbstractVector, b::AbstractVector) = evaluate(Euclidean(), a, b)
#evaluate{T <: Number}(dist::Euclidean, a::T, b::T) = abs(a-b)
#euclidean{T <: Number}(a::T, b::T) = evaluate(Euclidean(), a, b)
#

