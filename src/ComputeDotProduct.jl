#=
 @file ComputeDotProduct_ref.cpp

 HPCG routine
=#
using MPI
include("mytimer.jl")

#=
  Routine to compute the dot product of two vectors where:

  This is the reference dot-product implementation.  It _CANNOT_ be modified for the
  purposes of this benchmark.

  @param[in] n the number of vector elements (on this processor)
  @param[in] x, y the input vectors
  @param[in] result a pointer to scalar value, on exit will contain result.
  @param[out] time_allreduce the time it took to perform the communication between processes

  @return returns 0 upon success and non-zero otherwise

  @see ComputeDotProduct
=#
function ComputeDotProduct_ref(n, x, y, result, time_allreduce) 
  @assert(length(x)>=n) # Test vector lengths
  @assert(length(y)>=n)

  local_result = 0.0
  xv = x
  yv = y
  if yv==xv
    for i=1:n
	local_result += xv[i]*xv[i]
    end
    else 
    for i=1:n
	local_result += xv[i]*yv[i]
    end
  end

  #Use MPI's reduce function to collect all partial sums
  t0 = mytimer()
  global_result = 0.0
  MPI.Allreduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, MPI.COMM_WORLD)
  result = global_result
  time_allreduce += mytimer() - t0
  time_allreduce += 0.0
  result = local_result

  return 0
end
