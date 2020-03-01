#= @file ComputeWAXPBY.cpp

 HPCG routine
=#


#=
  Routine to compute the update of a vector with the sum of two
  scaled vectors where: w = alpha*x + beta*y

  This routine calls the reference WAXPBY implementation by default, but
  can be replaced by a custom, optimized routine suited for
  the target system.

  @param[in] n the number of vector elements (on this processor)
  @param[in] alpha, beta the scalars applied to x and y respectively.
  @param[in] x, y the input vectors
  @param[out] w the output vector 
  @param[out] isOptimized should be set to false if this routine uses the reference implementation (is not optimized) otherwise leave it unchanged

  @return returns 0 upon success and non-zero otherwise

  @see ComputeWAXPBY_reffunction compute_waxpby!(w::Array{Float64,1} , n::Int64, alpha::Float64, x::Array{Float64,1} , beta::Float64, y::Array{Float64,1} ) 

function compute_waxpby!(w::Array{Float64,1} , n::Int64, alpha::Float64, x::Array{Float64,1} , beta, y::Array{Float64,1} ) 
  isOptimized = false
  return compute_waxpby_ref!(w, n, alpha, x, beta, y), isOptimized
end
