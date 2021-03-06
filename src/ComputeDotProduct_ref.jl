#=
 @file ComputeDotProduct_ref.cpp

 HPCG routine
=#

#=
  Routine to compute the dot product of two vectors where:

  This is the reference dot-product implementation.  It _CANNOT_ be modified for the
  purposes of this benchmark.

  @param[in] n the number of vector elements (on this processor)
  @param[in] x, y the input vectors
  @param[in] result a pointer to scalar value, on exit will contain result.
  @param[out] time_allreduce the time it took to perform the communication between processes

  @return returns false upon success and true otherwise

  @see compute_dot_product
=#
function compute_dot_product_ref!(n::Int64, x::Array{Float64,1}, y::Array{Float64,1}, result::Float64, time_allreduce::Float64) 
  @assert(length(x)>=n) # Test vector lengths
  @assert(length(y)>=n)
  global_result = [0.0]
  local_result::Float64  = 0.0
#=  if y==x
     for i = 1:n
        local_result += x[i]*x[i]
    end
  else 
  =#
    for i = 1:n
      @fastmath @inbounds  local_result::Float64 += x[i]*y[i]
  	end
# end
  if MPI.Initialized() == true
 	#Use MPI's reduce function to collect all partial sums
  	t0::Float64 = time_ns()
#	global_result = [0.0]
  
	MPI.Allreduce!([local_result], global_result, 1, MPI.SUM, MPI.COMM_WORLD)
  	result = global_result[1]
  	time_allreduce::Float64 += time_ns() - t0
  return result, time_allreduce, 0
  else 
  	time_allreduce::Float64 += 0.0
   	result::Float64 = local_result
  return result, time_allreduce, 0
  end 
end
