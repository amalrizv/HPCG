#=
 @file ComputeSYMGS_ref.cpp

 HPCG routine
=#
include("ComputeMG_ref.jl")
include("ComputeSYMGS_ref.jl")
include("ComputeSPMV_ref.jl")
include("ComputeRestriction_ref.jl")
include("ComputeProlongation_ref.jl")

#=

  @param[in] A the known system matrix
  @param[in] r the input vector
  @param[inout] x On exit contains the result of the multigrid V-cycle with r as the RHS, x is the approximation to Ax = r.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeMG
=#
function ComputeMG_ref(A, r, x) 
  @assert(x.localLength==A.localNumberOfColumns) Make sure x contain space for halo values

  ZeroVector(x) #initialize x to zero

  ierr = 0
  if A.mgData!=0 #  Go to next coarse level if defined
    numberOfPresmootherSteps = A.mgData->numberOfPresmootherSteps
    for i=1: numberOfPresmootherSteps
	 ierr += ComputeSYMGS_ref(A, r, x)
    end
    if ierr!=0
	 return ierr
    end
    ierr = ComputeSPMV_ref(A, x, *A.mgData->Axf) 
    if ierr!=0
	 return ierr
    end
    # Perform restriction operation using simple injection
    ierr = ComputeRestriction_ref(A, r)  
    if ierr!=0 
	return ierr
    end
    ierr = ComputeMG_ref(*A.Ac,*A.mgData->rc, *A.mgData->xc)  
    if ierr!=0 
	return ierr
    end
    ierr = ComputeProlongation_ref(A, x)  
    if (ierr!=0) 
	return ierr
    end
    numberOfPostsmootherSteps = A.mgData->numberOfPostsmootherSteps
    for i= 1: numberOfPostsmootherSteps
	ierr += ComputeSYMGS_ref(A, r, x)
    end
    if ierr!=0 
	return ierr
    else 
    	ierr = ComputeSYMGS_ref(A, r, x)
    end 
    if ierr!=0 
	return ierr
    end
  end
  return 0
end

