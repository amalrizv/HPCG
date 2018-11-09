type  counter{T}
	length::T # number of prime factor counts(cannot exceed 32 or a 32-bit integer)
	 max_counts #maximum value for prime factor counts
	cur_counts #current prime factor counts
end
l = 33
a = Array{Int64,1}(33)
new = counter(l,a,a)
println(new)
