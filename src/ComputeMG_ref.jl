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
  if A.mgData == 0  #  Go to next coarse level if defined
      numberOfPresmootherSteps = A.mgData.numberOfPresmootherSteps

      # In second iteration Compute MG sends A.mgData.xc 
      # as inout and A.Ac and A.mgData.rc as inputs
      # For the second rank (rc and Ac inputs checked)
      # something is wrong with computation when rank is 1 
      for i = 1:numberOfPresmootherSteps
          ierr = compute_symgs_ref!(x, A, r )
      end
      if ierr!=0
          return ierr
      end
      ierr = compute_spmv_ref!(A.mgData.Axf, A, x) 
#      @show A.mgData.Axf[length(A.mgData.Axf)]
      #second iteration Axf is wrong for rank 1 
	
      if ierr!=0
          return ierr
      end

      # Perform restriction operation using simple injection
      ierr = compute_restriction_ref!(A, r)  
      if ierr!=0 
          return ierr
      end
      ierr, A.mgData.xc = compute_mg!(A.mgData.xc, A.Ac, A.mgData.rc)  
##############loops back one time################
      if ierr!=0 
          return ierr
      end

      ierr = compute_prolongation_ref!(x,A)  
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
      ierr = compute_symgs_ref!(x, A, r)
      if ierr!=0 
          return ierr
      end
  end

  return 0
end

