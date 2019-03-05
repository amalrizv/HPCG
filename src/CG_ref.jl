#=
 @file CG_ref.cpp
 HPCG routine
=#



include("hpcg.jl")

include("mytimer.jl")
include("ComputeSPMV_ref.jl")
include("ComputeMG_ref.jl")
include("ComputeDotProduct_ref.jl")
include("ComputeWAXPBY_ref.jl")



#=
  Reference routine to compute an approximate solution to Ax = b
  @param[inout] A    The known system matrix
  @param[inout] data The data structure with all necessary CG vectors preallocated
  @param[in]    b    The known right hand side vector
  @param[inout] x    On entry: the initial guess on exit: the new approximate solution
  @param[in]    max_iter  The maximum number of iterations to perform, even if tolerance is not met.
  @param[in]    tolerance The stopping criterion to assert convergence: if norm of residual is <= to tolerance.
  @param[out]   niters    The number of iterations actually performed.
  @param[out]   normr     The 2-norm of the residual vector after the last iteration.
  @param[out]   normr0    The 2-norm of the residual vector before the first iteration.
  @param[out]   times     The 7-element vector of the timing information accumulated during all of the iterations.
  @param[in]    doPreconditioning The flag to indicate whether the preconditioner should be invoked at each iteration.
  @return Returns zero on success and a non-zero value otherwise.
  @see CG()
=#
function CG_ref(const A, data, const b, x, const max_iter, const tolerance, niters, normr, normr0, times, doPreconditioning) 

  t_begin = time_ns()  # Start timing right away
  normr = 0.0
  double rtz = 0.0
  oldrtz = 0.0
  alpha = 0.0 
  beta = 0.0
  pAp = 0.0


  t0 = 0.0
  t1 = 0.0 
  t2 = 0.0 
  t3 = 0.0 
  t4 = 0.0
  t5 = 0.0
  t6 = 0.0

  nrow = A.localNumberOfRows

  r = data.r # Residual vector
  z = data.z # Preconditioned residual vector
  p = data.p # Direction vector (in MPI mode ncol>=nrow)
  Ap = data.Ap

  if !doPreconditioning && A.geom.rank==0
	@debug("WARNING: PERFORMING UNPRECONDITIONED ITERATIONS")
  end

#ifdef HPCG_DEBUG
  print_freq = 1
  if print_freq>50 
	print_freq=50
  end
  if print_freq<1
	print_freq=1
  end
#endif
  # p is of length ncols, copy x to p for sparse MV operation
  CopyVector(x, p)
  tic() 
	ComputeSPMV_ref(A, p, Ap)  
  t3 = toc() # Ap = A*p
  tic() 
  ComputeWAXPBY_ref(nrow, 1.0, b, -1.0, Ap, r) 
  t2 = toc() # r = b - Ax (x stored in p)
  tic() 
  ComputeDotProduct_ref(nrow, r, r, normr, t4)  
  t1 = toc()
  normr = sqrt(normr)
#ifdef HPCG_DEBUG
  if A.geom.rank==0 
     @debug("Initial Residual = $normr") 
  end
#endif

  # Record initial residual for convergence testing
  normr0 = normr

  # Start iterations
  while nomr/nomr0 > tolerance
    for (k=1:max_iter 
    	tic()
    	if doPreconditioning
      		ComputeMG_ref(A, r, z) # Apply preconditioner
    	else
      		ComputeWAXPBY_ref(nrow, 1.0, r, 0.0, r, z) # copy r to z (no preconditioning)
    	toc(t5) # Preconditioner apply time

    	if k == 1
      		CopyVector(z, p) toc(t2) # Copy Mr to p
      		tic() 
      		ComputeDotProduct_ref(nrow, r, z, rtz, t4) 
      		t1 = toc() # rtz = r'*z
    	else 
      		oldrtz = rtz
      		tic() 
      		ComputeDotProduct_ref(nrow, r, z, rtz, t4) 
      		t1 = toc() # rtz = r'*z
      		beta = rtz/oldrtz
      		tic() 
      		ComputeWAXPBY_ref(nrow, 1.0, z, beta, p, p)  
      		t2 = toc() # p = beta*p + z
    	end

    	tic() 
    	ComputeSPMV_ref(A, p, Ap) 
    	t3 = toc() # Ap = A*p
    	tic() 
    	ComputeDotProduct_ref(nrow, p, Ap, pAp, t4) 
    	t1 = toc() # alpha = p'*Ap
    	alpha = rtz/pAp
    	tic() 
    	ComputeWAXPBY_ref(nrow, 1.0, x, alpha, p, x)# x = x + alpha*p
    	ComputeWAXPBY_ref(nrow, 1.0, r, -alpha, Ap, r)  
    	t2 = toc()# r = r - alpha*Ap
    	tic() 
    	ComputeDotProduct_ref(nrow, r, r, normr, t4) 
    	t1 = toc()
    	normr = sqrt(normr)
#ifdef HPCG_DEBUG
    	if A.geom.rank==0 & k%print_freq == 0 || k == max_iter
      		@debug("Iteration = $k Scaled Residual = $(normr/normr0)")
    	end
#endif
   	niters = k
  end
 end

  # Store times
  times[1] += t1 # dot product time
  times[2] += t2 # WAXPBY time
  times[3] += t3 # SPMV time
  times[4] += t4 # AllReduce time
  times[5] += t5 # preconditioner apply time
##ifndef HPCG_NO_MPI
#  times[6] += t6 # exchange halo time
##endif
  times[0] += mytimer() - t_begin  # Total time. All done...
  return 0
end

