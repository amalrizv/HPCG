 #=@file TestNorms.hpp

 HPCG data structure
 =#


mutable struct TestNormsData
	values  #!< sample values
	mean    #!< mean of all sampes
	variance  #!< variance of mean
	samples   #!< number of samples
	pass      #!< pass/fail indicator
end


