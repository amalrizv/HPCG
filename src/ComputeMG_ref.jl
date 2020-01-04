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

function compute_mg_ref!(x, A, r) #sp_coarse passed 
  @assert(length(x)==A.localNumberOfColumns) #Make sure x contain space for halo values
  zero_fill!(x)
  #initialize x to zero

  if A.geom.rank == 0
	  fs =  open("mpi_ircv_0.txt", "a")
	  fp =  open("mpi_send_0.txt", "a")
	  sfs =  open("mg_symgs_call_order_0.txt", "a")
  else
 	  fs =  open("mpi_ircv_1.txt", "a")
	  fp =  open("mpi_send_1.txt", "a")
	  sfs =  open("mg_symgs_call_order_1.txt", "a")
  end

  if A.mgData.init == true
	  println(sfs, "A.mgData != 0; PreCoarseningSmoothening ")

      numberOfPresmootherSteps = A.mgData.numberOfPresmootherSteps
      for i = 1:numberOfPresmootherSteps

		# println(fs, "symgs x Called from MG")
	 	# println(fp, "symgs x Called from MG")
      	ierr = compute_symgs_ref!(x, A, r )
      	if ierr!=0
          return ierr
      	end
      end

	  #println(fs, "spmv A.mgData.Axf Called from MG")
	  #println(fp, "spmv A.mgData.Axf Called from MG")

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

	  #println(fs, "Recursive MG A.mgData.xc from MG")
	  #println(fp, "Recursive MG A.mgData.xc from MG")
	  #A.mgData.xc is wrong input for SYMGS after recursive MG
	  println(sfs, "A.mgData != 0 RECURSIVE CALL TO mg")
	  close(sfs)
      ierr = compute_mg_ref!(A.mgData.xc, A.Ac, A.mgData.rc)  
      if ierr!=0 
          return ierr
      end

      ierr = compute_prolongation_ref!(x,A)  
      if (ierr!=0) 
          return ierr
      end

      numberOfPostsmootherSteps = A.mgData.numberOfPostsmootherSteps
	if A.geom.rank == 0
	 	sfs =  open("mg_symgs_call_order_0.txt", "a")
  	else
		sfs =  open("mg_symgs_call_order_1.txt", "a")	
	end

	  println(sfs, "A.mgData != 0; PostCoarseningSmoothening ")
      for i= 1: numberOfPostsmootherSteps

	#	    println(fs, "symgs x [seclast symgs call of mg] Called from MG")
	#		println(fp, "symgs x [seclast symgs call of mg] Called from MG")
    		ierr = compute_symgs_ref!(x,A, r)
      	if ierr!=0 
        	  return ierr
      	end

      end
  else 
	  println(sfs, "A.mgData = 0 ")

	 # println(fs, "symgs x [last symgs call of mg] Called from MG")
	 # println(fp, "symgs x [last symgs call of mg] Called from MG")
      ierr = compute_symgs_ref!(x, A, r)
      if ierr!=0 
          return ierr
      end

  end
  close(fs)
  close(fp)
  close(sfs)
  
  return ierr
end

