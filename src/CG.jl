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
function cg!(A, data, b, x, max_iter, tolerance, times, doPreconditioning) 

  t_begin = time_ns()  # Start timing right away
  normr = 0.0
  normr0 = 0.0
  niters = 0
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


  p[1:length(x)] = x
  t3t 		= time_ns() 
  ierr = compute_spmv!(Ap, A, p) 
  t3 		= time_ns()-t3t 

  #Ap 	= A*p

  t2t 	= time_ns()
   ierr, A.is_waxpby_optimized  = compute_waxpby!(r, nrow, 1.0, b, -1.0, Ap)

  t2 	= time_ns()-t2t 

  # r = b - Ax (x stored in p)

  t1t 	= time_ns()
  normr, t4, ierr = compute_dot_product!(nrow, r, r)
  t1 	= time_ns()-t1t

  normr = sqrt(normr)
#ifdef HPCG_DEBUG
  if A.geom.rank == 0 
  	@debug("Initial Residual = $normr")
  end
#endif

  # Record initial residual for convergence testing
  normr0 = normr


  # Start iterations
  for k	= 1:max_iter
	  if  normr/normr0 > tolerance && iszero(normr) == false
		       		t5t = time_ns()
    		if doPreconditioning == true
    #symgs->exchnagehalo
      			ierr = compute_mg!(z,A, r) # Apply preconditioner
    		else
      			z[1:length(r)] = r # copy r to z (no preconditioning)
    		end
    		t5 	= time_ns()- t5t # Preconditioner apply time

    		if k == 1
      			t2t   	= time_ns()
      			ierr, A.is_waxpby_optimized  =  compute_waxpby!(p, nrow, 1.0, z, 0.0, z)
      			t2 	= t2+time_ns()-t2t # Copy Mr to p
      			t1t 	= time_ns()
      			rtz, t4, ierr  = compute_dot_product!(nrow, r, z)
      			t1 	= t1+time_ns()- t1t # rtz = r'*z
   			else 
      			oldrtz 	= rtz
      			t1t 	= time_ns()
      		    rtz, t4, ierr  = compute_dot_product!( nrow, r, z) 
      			t1 	= t1+time_ns()-t1t # rtz = r'*z
      			beta 	= rtz/oldrtz
      			t2t 	= time_ns()
      			ierr, A.is_waxpby_optimized 	=compute_waxpby!(p, nrow, 1.0, z, beta, p)  
      			t2 	= time_ns()-t2t+t2 # p = beta*p + z
  	 		end
			
    		t3t 	= time_ns()
    		ierr = compute_spmv!(Ap, A, p) 
    		t3	= t3+time_ns()- t3t # Ap = A*p

    		t1t 	= time_ns()
    		pAp, t4, ierr  = compute_dot_product!(nrow, p, Ap) 
    		t1 	= time_ns()-t1t+t1 # alpha = p'*Ap
  		
			if pAp == 0
				alpha = 1
			else
	    		alpha 	= rtz/pAp
			end
    		t2t  	= time_ns() 

 
		
    		ierr, A.is_waxpby_optimized, = compute_waxpby!(x, nrow, 1.0, x, alpha, p)# x = x + alpha*p
  
    		ierr, A.is_waxpby_optimized  = compute_waxpby!(r, nrow, 1.0, r, -alpha, Ap)
    		t2 	= time_ns()- t2t +t2# r = r - alpha*Ap
	
    		t1t 	= time_ns()
            normr, t4, ierr  = compute_dot_product!(nrow, r, r)
    		t1 = t1+time_ns()- t1t
	
    		normr = sqrt(normr)

    		if A.geom.rank==0 && (k%print_freq == 0 || k == max_iter)
				@debug("Iteration = $k    Scaled Residual = $(normr/normr0)")
			end
    		niters = k
			#@show normr
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
	  
 ierr=0
 	#= 
	data.r  = r # Residual vector
    data.z  = z # Preconditioned residual vector
    data.p  = p # Direction vector (in MPI mode ncol>=nrow)
    data.Ap = Ap
	=#
	#@show "About to exit CG"
	#@show normr 
	A.is_dot_prod_optimized = false
 return niters, normr, normr0, ierr
end
