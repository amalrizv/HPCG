#=
 @file ComputeSYMGS_ref.cpp

 HPCG routine
=#
include("ComputeSPMV_ref.jl")
include("ComputeRestriction_ref.jl")
include("ComputeSYMGS_ref.jl")
include("ComputeProlongation_ref.jl")

#=

  @param[in] A the known system matrix
  @param[in] r the input vector
  @param[inout] x On exit contains the result of the multigrid V-cycle with r as the RHS, x is the approximation to Ax = r.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeMG
=#

function compute_mg_ref!(x, A, r) #sp_coarse passed 
  @assert(length(x)==A.localNumberOfColumns) #Make sure x contain space for halo values

  x = zeros(length(x)) #initialize x to zero

  ierr = 0

  # if mgdata is even defined this will result in a non zero value 
  if length(A.mgData.Axf) != 0 #  Go to next coarse level if defined
      numberOfPresmootherSteps = A.mgData.numberOfPresmootherSteps

      for i = 1:numberOfPresmootherSteps
          ierr, x = compute_symgs_ref!(x, A, r )
      end
      if ierr!=0
          return ierr
      end

      ierr, A.mgData.Axf = compute_spmv_ref!(A.mgData.Axf, A, x) 
      if ierr!=0
          return ierr
      end

      # Perform restriction operation using simple injection
      ierr, A= compute_restriction_ref!(A, r)  
      if ierr!=0 
          return ierr
      end
      ierr, A.mgData.xc = compute_mg!(A.mgData.xc, A.Ac, A.mgData.rc)  
      if ierr!=0 
          return ierr
      end

      ierr, x = compute_prolongation_ref!(x,A)  
      if (ierr!=0) 
          return ierr
      end

      numberOfPostsmootherSteps = A.mgData.numberOfPostsmootherSteps
      for i= 1: numberOfPostsmootherSteps
          ierr, x = compute_symgs_ref!(x,A, r)
      end

      if ierr!=0 
          return ierr
      end
  else 
      ierr, x = compute_symgs_ref!(x, A, r)
      if ierr!=0 
          return ierr
      end
  end

  return 0, x
end

