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

function compute_mg_ref!(A, r, x) #sp_coarse passed 

  #@assert(length(x)==A.localNumberOfColumns) #Make sure x contain space for halo values

  x = zeros(length(x)) #initialize x to zero

  ierr = 0

  println("rx")
  # if mgdata is even defined this will result in a non zero value 
  if A.mgData != 0 #  Go to next coarse level if defined
      numberOfPresmootherSteps = A.mgData.numberOfPresmootherSteps

      for i = 1:numberOfPresmootherSteps
          ierr += compute_symgs_ref(A, r, x)
      end
      if ierr!=0
          return ierr
      end

      @show length(A.mgData.Axf)
      ierr = compute_spmv_ref(A, x, A.mgData.Axf) 
      if ierr!=0
          return ierr
      end

      # Perform restriction operation using simple injection
      ierr = compute_restriction_ref(A, r)  
      if ierr!=0 
          return ierr
      end

      ierr = compute_mg!(A.Ac, A.mgData.rc, A.mgData.xc)  
      if ierr!=0 
          return ierr
      end

      ierr = compute_prolongation_ref!(A, x)  
      if (ierr!=0) 
          return ierr
      end

      numberOfPostsmootherSteps = A.mgData.numberOfPostsmootherSteps
      for i= 1: numberOfPostsmootherSteps
          ierr += compute_symgs_ref(A, r, x)
      end

      if ierr!=0 
          return ierr
      end
  else 
      ierr = compute_symgs_ref(A, r, x)
      if ierr!=0 
          return ierr
      end
  end

  return 0
end

