#=
 @file TestNorms.cpp

 HPCG routine
=#

include("TestNorms_struct.jl")

#=
  Computes the mean and standard deviation of the array of norm results.

  @param[in] testnorms_data data structure with the results of norm test

  @return Returns 0 upon success or non-zero otherwise
=#
function test_norms(testnorms_data) 

 mean_delta = 0.0

 for i= 1:testnorms_data.samples
	mean_delta += (testnorms_data.values[i] - testnorms_data.values[1])
 end

 mean = testnorms_data.values[1] + mean_delta/testnorms_data.samples
 testnorms_data.mean = mean

 #Compute variance
 sumdiff = 0.0

 for i= 1:testnorms_data.samples 
	sumdiff += (testnorms_data.values[i] - mean) * (testnorms_data.values[i] - mean)
 end

 testnorms_data.variance = sumdiff/testnorms_data.samples

 #Determine if variation is sufficiently small to declare success
 testnorms_data.pass = (testnorms_data.variance<1.0e-6)

 return 0
end 
