#=
 @file ComputeDotProduct_ref.cpp

 HPCG routine
=#

using MPI
#include("mytimer.jl")

#=
  Routine to compute the dot product of two vectors where:

  This is the reference dot-product implementation.  It _CANNOT_ be modified for the
  purposes of this benchmark.

  @param[in] n the number of vector elements (on this processor)
  @param[in] x, y the input vectors
  @param[inout] result a pointer to scalar value, on exit will contain result.
  @param[out] time_allreduce the time it took to perform the communication between processes

  @return returns false upon success and true otherwise

  @see compute_dot_product
=#
function compute_dot_product_ref!(result, time_allreduce, n, x, y) 

  @assert(length(x)>=n) # Test vector lengths
  @assert(length(y)>=n)

  local_result = 0.0

  if x == y
    for i=1:n
        local_result += x[i]*x[i]
    end
  else 
    for i=1:n
        local_result += x[i]*y[i]
    end
  end

@static if MPI.Initialized()
  #Use MPI's reduce function to collect all partial sums
  t0 = time_ns()
  global_result = 0.0
  MPI.Allreduce!(local_result, global_result, 1, MPI.DOUBLE, MPI.SUM, MPI.COMM_WORLD)
  result = global_result
  time_allreduce += time_ns() - t0
else 
  time_allreduce += 0.0
  result = local_result
end

  return false

end