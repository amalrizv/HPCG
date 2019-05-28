#=
@file ComputeResidual.cpp

 HPCG routine
=#
using MPI

include("hpcg.jl")


#=
  Routine to compute the inf-norm difference between two vectors where:

  @param[in]  n        number of vector elements (local to this processor)
  @param[in]  v1, v2   input vectors
  @param[out] residual pointer to scalar value on exit, will contain result: inf-norm difference

  @return Returns zero on success and a non-zero value otherwise.
=#
function compute_residual(n, v1, v2, residual) 

  v1v = v1
  v2v = v2
  local_residual = 0.0

     threadlocal_residual = 0.0
    for i=1:n
       diff = (v1v[i] - v2v[i])
      if diff > threadlocal_residual 
	threadlocal_residual = diff
      end
   
    
      if threadlocal_residual>local_residual
	 local_residual = threadlocal_residual
      end
    
  end
  for i=1:n
     diff = (v1v[i] - v2v[i])
    if diff > local_residual
	 local_residual = diff
    end
    @debug(" Computed, exact, diff = $v1v[i] $v2v[i] $diff")
   end

  if MPI.Initialized()== true
	  # Use MPI's reduce function to collect all partial sums
  	global_residual = 0
  	global_residual = MPI.Allreduce(local_residual, MPI.MAX, MPI.COMM_WORLD)
  	residual = global_residual
  else
 	 residual = local_residual
  end
  return 0
end
