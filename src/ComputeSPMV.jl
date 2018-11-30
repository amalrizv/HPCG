
#=
 @file ComputeSPMV.cpp

 HPCG routine
=#

include("ComputeSPMV_ref.jl")

#=
  Routine to compute sparse matrix vector product y = Ax where:
  Precondition: First call exchange_externals to get off-processor values of x

  This routine calls the reference SpMV implementation by default, but
  can be replaced by a custom, optimized routine suited for
  the target system.

  @param[in]  A the known system matrix
  @param[in]  x the known vector
  @param[out] y the On exit contains the result: Ax.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeSPMV_ref
=#
function ComputeSPMV(A, x,y) 

  # This line and the next two lines should be removed and your version of ComputeSPMV should be used.
  A.isSpmvOptimized = false
  return ComputeSPMV_ref(A, x, y)
end
