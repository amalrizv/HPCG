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
function cg!(A, data, b, x, max_iter, tolerance, niters, normr, normr0, times, doPreconditioning) 
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
  x = p
  t3t = time_ns() 
  flag, Ap = compute_spmv!(Ap, A, p) 
  t3 = time_ns()-t3t 
  #Ap = A*p
  t2t = time_ns()
  flag, r = compute_waxpby!(r, nrow, 1.0, b, -1.0, Ap, A.is_waxpby_optimized)
  t2 = time_ns()-t2t 
  # r = b - Ax (x stored in p)
  t1t = time_ns()
  flag, normr = compute_dot_product!(normr, t4, nrow, r, r, A.is_dot_prod_optimized)
  t1 = time_ns()-t1t
  normr = sqrt(normr)
  @show normr
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
    		t5t = time_ns()
    		if doPreconditioning
      			flag, z = compute_mg!(z,A, r) # Apply preconditioner
    		else
      			z = r # copy r to z (no preconditioning)
    		end
    		t5 = time_ns()- t5t # Preconditioner apply time

    		if k == 1
      			t2t = time_ns()
      			flag, p = compute_waxpby!(p, nrow, 1.0, z, 0.0, z, A.is_waxpby_optimized)
      			t2 = t2+time_ns()-t2t # Copy Mr to p
      			t1t = time_ns()
      			flag, rtz = compute_dot_product!(rtz, t4, nrow, r, z, A.is_dot_prod_optimized)
      			t1 = t1+time_ns()- t1t # rtz = r'*z
   		else 
      			oldrtz = rtz
      			t1t = time_ns()
      			flag, rtz = compute_dot_product!(rtz, t4, nrow, r, z, A.is_dot_prod_optimized) 
      			t1 = t1+time_ns()-t1t # rtz = r'*z
      			beta = rtz/oldrtz
      			t2t = time_ns()
      			flag, p = compute_waxpby!(p, nrow, 1.0, z, beta, p, A.is_waxpby_optimized)  
      			t2 = time_ns()-t2t+t2 # p = beta*p + z
   		end

    		t3t = time_ns()
    		flag, Ap = compute_spmv!(Ap, A, p) 
    		t3 = t3+time_ns()- t3t # Ap = A*p
    		t1t = time_ns()
    		flag, pAp = compute_dot_product!(pAp, t4, nrow, p, Ap, A.is_dot_prod_optimized) 
    		t1 = time_ns()-t1t+t1 # alpha = p'*Ap
    		alpha = rtz/pAp
    		t2t  = time_ns() 
    		flag, x = compute_waxpby!(x, nrow, 1.0, x, alpha, p, A.is_waxpby_optimized)# x = x + alpha*p
    		flag, r = compute_waxpby!(r, nrow, 1.0, r, -alpha, Ap, A.is_waxpby_optimized)
    		t2 = time_ns()- t2t +t2# r = r - alpha*Ap
    		t1t = time_ns()
                flag, normr = compute_dot_product!(normr, t4, nrow, r, r, A.is_dot_prod_optimized)
    		t1 = t1+time_ns()- t1t
    		normr = sqrt(normr)
#ifdef HPCG_DEBUG
    		if A.geom.rank==0 && (k%print_freq == 0 || k == max_iter)
      			@debug("Iteration = ",k,"   Scaled Residual = ",normr/normr0, "\n")
		end
#endif
    		niters = k
  	end
  end

  t0 = time_ns() - t_begin  # Total time. All done...
  times[1] += t0
  times[2] += t1
  times[3] += t2
  times[4] += t3
  times[5] += t4
  times[6] += t5
##ifndef HPCG_NO_MPI
## times[7] = t6
##else
##endif
  
  return 0 ,A , data , x , niters, normr , normr0 , times
end
