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

function ComputeMG(AAAA, r, x, sp) #sp_coarse passed 
  if sp == 0
   AAA = AAAA.Ac
  else
   AAA = AAAA.sp_matrix 		# sp_anx structure
  end
  AA = AAA.sp_matrix		#sp_matrix structure
  #@assert(length(x)==AAA.localNumberOfCols) #Make sure x contain space for halo values

  x = zeros(length(x)) #initialize x to zero

  ierr = 0
  # if mgdata is even defined this will result in a non zero value 
  if sp == 1 #  Go to next coarse level if defined
    numberOfPresmootherSteps = AAAA.mgData.numberOfPresmootherSteps
    for i=1: numberOfPresmootherSteps
	 ierr += ComputeSYMGS_ref(AAA, r, x)
    end
    if ierr!=0
	 return ierr
    end
    ierr = ComputeSPMV_ref(AAA, x, AAAA.mgData.Axf) 
    if ierr!=0
	 return ierr
    end
    # Perform restriction operation using simple injection
    ierr = ComputeRestriction_ref(AAAA, r)  
    if ierr!=0 
	return ierr
    end
    ierr = ComputeMG(AAAA,AAAA.mgData.rc, AAAA.mgData.xc, 0)  
    if ierr!=0 
	return ierr
    end
    ierr = ComputeProlongation_ref(AAAA, x)  
    if (ierr!=0) 
	return ierr
    end
    numberOfPostsmootherSteps = AAAA.mgData.numberOfPostsmootherSteps
    for i= 1: numberOfPostsmootherSteps
	ierr += ComputeSYMGS_ref(AAA, r, x)
    end
    if ierr!=0 
	return ierr
    end
  else 
    ierr = ComputeSYMGS_ref(AAA, r, x)
    if ierr!=0 
	return ierr
    end
  end
  return 0
end

