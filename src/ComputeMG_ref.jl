#=
 @file ComputeSYMGS_ref.cpp

 HPCG routine
=#
include("appendx.jl")
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

function compute_mg_ref!(x, A, r, ierr) #sp_coarse passed 
  @assert(length(x)==A.localNumberOfColumns) #Make sure x contain space for halo values
  zero_fill!(x)
  #initialize x to zero


  # if mgdata is even defined this will result in a non zero value 
  if A.mgData == 0  #  Go to next coarse level if defined
      numberOfPresmootherSteps = A.mgData.numberOfPresmootherSteps

      for i = 1:numberOfPresmootherSteps
          compute_symgs_ref!(x, A, r , ierr)
      end
      if ierr!=0
          return ierr
      end

      ierr = compute_spmv_ref!(A.mgData.Axf, A, x) 
      if ierr!=0
          return ierr
      end

      # Perform restriction operation using simple injection
      compute_restriction_ref!(A, r, ierr)  
      if ierr!=0 
          return ierr
      end

      compute_mg!(A.mgData.xc, A.Ac, A.mgData.rc)  
      if ierr!=0 
          return ierr
      end

      compute_prolongation_ref!(x,A, ierr)  
      if (ierr!=0) 
          return ierr
      end

      numberOfPostsmootherSteps = A.mgData.numberOfPostsmootherSteps

      for i= 1: numberOfPostsmootherSteps
    		compute_symgs_ref!(x,A, r, ierr)
      end
      if ierr!=0 
          return ierr
      end

  else 

      compute_symgs_ref!(x, A, r, ierr)
      if ierr!=0 
          return ierr
      end

  end
  
  ierr = 0
end

