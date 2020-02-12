#=
 @file ComputeSYMGS_ref.cpp

 HPCG routine
=#


#=

  @param[in] A the known system matrix
  @param[in] r the input vector
  @param[inout] x On exit contains the result of the multigrid V-cycle with r as the RHS, x is the approximation to Ax = r.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeMG
=#

function compute_mg_ref!(x, A, r) #sp_coarse passed 
  @assert(length(x)==A.localNumberOfColumns) #Make sure x contain space for halo values
  zero_fill!(x)
  #initialize x to zero

  
  if A.mgData.init == true

      numberOfPresmootherSteps = A.mgData.numberOfPresmootherSteps
      for i = 1:numberOfPresmootherSteps

      	ierr = compute_symgs_ref!(x, A, r )
      	if ierr!=0
          return ierr
      	end
      end


	  # x vector in the first call of line 56 should be 0
	  # instead it is NaN. This is the first instance 
	  # of ExchangeHalo discrepancy
	  
      ierr = compute_spmv_ref!(A.mgData.Axf, A, x) 
      if ierr!=0
          return ierr
      end

      # Perform restriction operation using simple injection
      ierr = compute_restriction_ref!(A, r)  
      if ierr!=0 
          return ierr
      end

	  #A.mgData.xc is wrong input for SYMGS after recursive MG
      ierr = compute_mg_ref!(A.mgData.xc, A.Ac, A.mgData.rc)  
      if ierr!=0 
          return ierr
      end

      ierr = compute_prolongation_ref!(x,A)  
      if (ierr!=0) 
          return ierr
      end

      numberOfPostsmootherSteps = A.mgData.numberOfPostsmootherSteps

      for i= 1: numberOfPostsmootherSteps

    		ierr = compute_symgs_ref!(x,A, r)
      	if ierr!=0 
        	  return ierr
      	end

      end
  else 

      ierr = compute_symgs_ref!(x, A, r)
      if ierr!=0 
          return ierr
      end

  end
   return ierr
end

