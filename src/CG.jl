#@HEADER
# ***************************************************
#
# HPCG: High Performance Conjugate Gradient Benchmark
#
# Contact:
# Michael A. Heroux ( maherou@sandia.gov)
# Jack Dongarra     (dongarra@eecs.utk.edu)
# Piotr Luszczek    (luszczek@eecs.utk.edu)
#
# ***************************************************
#@HEADER

#=
 @file CG.cpp
 HPCG routine
=#

#include <fstream>

#include <cmath>

include("hpcg.jl")


#include("mytimer.jl")
include("ComputeSPMV.jl")
include("ComputeMG.jl")
include("ComputeDotProduct.jl")
include("ComputeWAXPBY.jl")



#=
  Routine to compute an approximate solution to Ax = b
  @param[in]    geom The description of the problem's geometry.
  @param[inout] A    The known system matrix
  @param[inout] data The data structure with all necessary CG vectors preallocated
  @param[in]    b    The known right hand side vector
  @param[inout] x    On entry: the initial guess; on exit: the new approximate solution
  @param[in]    max_iter  The maximum number of iterations to perform, even if tolerance is not met.
  @param[in]    tolerance The stopping criterion to assert convergence: if norm of residual is <= to tolerance.
  @param[out]   niters    The number of iterations actually performed.
  @param[out]   normr     The 2-norm of the residual vector after the last iteration.
  @param[out]   normr0    The 2-norm of the residual vector before the first iteration.
  @param[out]   times     The 7-element vector of the timing information accumulated during all of the iterations.
  @param[in]    doPreconditioning The flag to indicate whether the preconditioner should be invoked at each iteration.
  @return Returns zero on success and a non-zero value otherwise.
  @see CG_ref()
=#
function  CG(A, data, b, x, max_iter, tolerance, niters, normr, normr0, times, doPreconditioning) 

  t_begin = time_ns()  # Start timing right away
  normr = 0.0
  rtz = 0.0
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
##ifndef HPCG_NO_MPI
# t6 = 0.0
##endif
  nrow = A.localNumberOfRow 
  r = data.r # Residual vector
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
  x = p
  tic() 
  ComputeSPMV(A, p, Ap) 
  t3 = toc() 
  #Ap = A*p
  tic() 
  ComputeWAXPBY(nrow, 1.0, b, -1.0, Ap, r, A.isWaxpbyOptimized)
  t2 = toc(t2) 
  # r = b - Ax (x stored in p)
  tic() 
  ComputeDotProduct(nrow, r, r, normr, t4, A.isDotProductOptimized)
  t1 = toc()
  normr = sqrt(normr)
#ifdef HPCG_DEBUG
  if A.geom.rank==0 
  	@debug("Initial Residual = ",normr,"\n")
  end
#endif

  # Record initial residual for convergence testing
  normr0 = normr

  # Start iterations
  while normr/normr0 > tolerance
  	for k=1:max_iter+1
    		tic()
    		if doPreconditioning
      			ComputeMG(A, r, z) # Apply preconditioner
    		else
      			z = r # copy r to z (no preconditioning)
    		end
    		toc(t5) # Preconditioner apply time

    		if k == 1
      			tic() 
      			ComputeWAXPBY(nrow, 1.0, z, 0.0, z, p, A.isWaxpbyOptimized)
      			t2 = toc() # Copy Mr to p
      			tic()
      			ComputeDotProduct(nrow, r, z, rtz, t4, A.isDotProductOptimized)
      			t1 = toc() # rtz = r'*z
   		else 
      			oldrtz = rtz
      			tic() 
      			ComputeDotProduct(nrow, r, z, rtz, t4, A.isDotProductOptimized) 
      			t1 = toc() # rtz = r'*z
      			beta = rtz/oldrtz
      			tic() 
      			ComputeWAXPBY(nrow, 1.0, z, beta, p, p, A.isWaxpbyOptimized)  
      			t2 = toc() # p = beta*p + z
   		end

    		tic() 
    		ComputeSPMV(A, p, Ap) 
    		t3 = toc() # Ap = A*p
    		tic() 
    		ComputeDotProduct(nrow, p, Ap, pAp, t4, A.isDotProductOptimized) 
    		t1 = toc() # alpha = p'*Ap
    		alpha = rtz/pAp
    		tic() 
    		ComputeWAXPBY(nrow, 1.0, x, alpha, p, x, A.isWaxpbyOptimized)# x = x + alpha*p
    		ComputeWAXPBY(nrow, 1.0, r, -alpha, Ap, r, A.isWaxpbyOptimized)
    		t2 = toc()# r = r - alpha*Ap
    		tic() 
		ComputeDotProduct(nrow, r, r, normr, t4, A.isDotProductOptimized)
    		t1 = toc()
    		normr = sqrt(normr)
#ifdef HPCG_DEBUG
    		if A.geom.rank==0 && (k%print_freq == 0 || k == max_iter)
      			@debug("Iteration = ",k,"   Scaled Residual = ",normr/normr0, "\n")
		end
#endif
    		niters = k
  	end
  end
  # Store times
  times[1] += t1 # dot-product time
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

