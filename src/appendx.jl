function zero_fill!(x::Array{Float64,1})
	x[:]= [0.0 for i in 1:length(x)]
end
