#=
@file ComputeResidual.cpp

 HPCG routine
=#
using MPI



#=
  Routine to compute the inf-norm difference between two vectors where:

  @param[in]  n        number of vector elements (local to this processor)
  @param[in]  v1, v2   input vectors
  @param[out] residual pointer to scalar value on exit, will contain result: inf-norm difference

  @return Returns zero on success and a non-zero value otherwise.
=#
function compute_residual!(n, v1, v2) 

  local_residual = 0.0
      
  for i=1:n # No threading
     diff = abs(v1[i] - v2[i])
   	 if diff > local_residual
	     local_residual = diff
     end
	 @debug(" Computed, exact, diff = $(v1[i]), $(v2[i]),  $(diff)")
  end

  if MPI.Initialized()== true
	  # Use MPI's reduce function to collect all partial sums
  	global_residual = 0.0
  	global_residual = MPI.Allreduce(local_residual, MPI.MAX, MPI.COMM_WORLD)
  	residual = global_residual
  else
 	 residual = local_residual
  end
  #@show local_residual
  #@show residual
  ierr = 0
  return residual, 0
end
