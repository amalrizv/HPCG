 #=@file TestNorms.hpp

 HPCG data structure
 =#


mutable struct TestNormsData
	values::Array{Float64}  #!< sample values
	mean::{Float64}    #!< mean of all sampes
	variance::Float64  #!< variance of mean
	samples::Int64   #!< number of samples
	pass::Bool      #!< pass/fail indicator
end


