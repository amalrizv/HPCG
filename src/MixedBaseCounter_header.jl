mutable struct MixedBaseCounter 
    length::Int64 # number of prime factor counts (cannot exceed 32 for a 32-bit integer)
    max_counts::Array{Int64}(undef, 32+1) #  maximum value for prime factor counts
    cur_counts::Array{Int64}(undef, 32+1) #  current prime factor counts
end
