module distances
using ArrayViews

export 
	Euclidean,
	euclidean,
	PreMetric,
	SemiMetric,
	Metric,
	result_type,
	evaluate

include("metrics.jl")

end


