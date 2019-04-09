#= @file ComputeWAXPBY.cpp

 HPCG routine
=#

include("ComputeWAXPBY_ref.jl")

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

  @see ComputeWAXPBY_ref
=#
function compute_waxpby!(w, n, alpha, x, beta, y, isOptimized) 

  # This line and the next two lines should be removed and your version of ComputeWAXPBY should be used.
  isOptimized = false
  return compute_waxpby_ref!(w, n, alpha, x, beta, y)

end
