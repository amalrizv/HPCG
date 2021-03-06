#=
 @file CG_ref.cpp
 HPCG routine
=#
#=
  Reference routine to compute an approximate solution to Ax = b
  @param[inout] A    The known system matrix
  @param[inout] data The data structure with all necessary CG vectors preallocated
  @param[in]    b    The known right hand side vector:
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

function cg_ref!(A::HPCGSparseMatrix , data::CGData , b::Array{Float64,1} , x::Array{Float64,1} , max_iter::Int64 , tolerance::Float64 , times::Array{Float64}, doPreconditioning::Bool) 

    t_begin = time_ns()  # Start timing right away
    normr   = 0.0
    rtz     = 0.0
    oldrtz  = 0.0
    alpha   = 0.0
    beta    = 0.0
    pAp     = 0.0
    niters = 0



    t0 = 0.0
    t1 = 0.0 
    t2 = 0.0 
    t3 = 0.0 
    t4 = 0.0
    t5 = 0.0
 

    #t6 = 0.0

    nrow = A.localNumberOfRows
    r  = data.r # Residual vector
    z  = data.z # Preconditioned residual vector
    p  = data.p # Direction vector (in MPI mode ncol>=nrow)
    Ap = data.Ap
 
    if !doPreconditioning && A.geom.rank==0
        @debug("WARNING: PERFORMING UNPRECONDITIONED ITERATIONS")
    end

    print_freq = 1

    if print_freq > 50 
        print_freq = 50
    end

    if print_freq<1
        print_freq = 1
    end

    # p is of length ncols, copy x to p for sparse MV operation
    x_len 	= length(x)
    p[1:x_len] 	= x

    t3t 	= time_ns()
    ierr = compute_spmv_ref!(Ap, A , p)  


    t3  += time_ns()-t3t # Ap = A*p
	t2t = time_ns()

    ierr =    compute_waxpby_ref!(r, nrow, 1.0, b, -1.0, Ap) 
    t2  += time_ns()- t2t # r = b - Ax (x stored in p)
    t1t = time_ns()
    #Profile.clear()
	normr, t4, ierr = compute_dot_product_ref!(nrow, r, r, normr, t4)
	#pprof(web=false)

    t1    += time_ns()- t1t
    normr = sqrt(normr)

    @debug("Initial Residual = $normr") 

    # Record initial residual for convergence testing
    normr0 = normr

    # Start iterations
    for k=1:max_iter 
    	if normr/normr0 > tolerance
            t5t = time_ns()
            if doPreconditioning == true
				#true passed from main 
                ierr = compute_mg_ref!(z,A, r ) # Apply preconditioner
	    else
                ierr = compute_waxpby_ref!(z, nrow, 1.0, r, 0.0, r) # copy r to z (no preconditioning)
            end
		
            t5 += time_ns()- t5t # Preconditioner apply time

            if k == 1

                p[1:length(z)] 	= z 
                t2	+= time_ns()-t5t # Copy Mr to p
                t1t 	= time_ns() 
		rtz, t4, ierr 	= compute_dot_product_ref!(nrow, r, z, rtz, t4) 
		t1 	+= time_ns()- t1t # rtz = r'*z

            else 

                oldrtz = rtz

				# rtz = r'*z
                t1t 	= time_ns() 
                rtz , t4, ierr 	= compute_dot_product_ref!(nrow, r, z, rtz, t4) 
                t1 	+= time_ns()- t1t 

                beta 	= rtz/oldrtz
				# p = beta*p + z
                t2t 	= time_ns() 
                ierr 	= compute_waxpby_ref!(p, nrow, 1.0, z, beta, p)  
                t2 	+= time_ns()- t2t 
            end


			# Ap = A*p
            t3t   	= time_ns() 
            ierr 	= compute_spmv_ref!(Ap,  A, p) 
            t3   	+= time_ns()- t3t 

			# alpha = p'*Ap
			t1t   			= time_ns()
            pAp, t4, ierr 	= compute_dot_product_ref!(nrow, p, Ap, pAp, t4) 
            t1    			+= time_ns()- t1t

			if pAp == 0 
				alpha = 1
			else
            	alpha 	= rtz/pAp
		    end
			
			# x = x + alpha*p ; r = r - alpha*Ap
            t2t   	= time_ns()
			ierr 	= compute_waxpby_ref!(x, nrow, 1.0, x, alpha, p)# x = x + alpha*p
			
            ierr 	= compute_waxpby_ref!(r, nrow, 1.0, r, -alpha, Ap)  
            t2    	+= time_ns()- t2t

            t1t   			=	 time_ns()
            normr, t4, ierr = compute_dot_product_ref!(nrow, r,r, normr, t4) 
            t1    			+= time_ns()- t1t

            normr 	= sqrt(normr)

            if A.geom.rank==0 & k%print_freq == 0 || k == max_iter
                @debug("Iteration = $k Scaled Residual = $(normr/normr0)")
            end

            niters = k
		else
			break
		end # end of tolerance check 
	end # end of for loop

    # Store times

    t0        = time_ns() - t_begin  # Total time. All done...
    times_add = Array{Any}(undef,length(times))

	times_add[1] = t0+times[1]
	times_add[2] = t1+times[2]
	times_add[3] = t2+times[3]
	times_add[4] = t3+times[4]
	times_add[5] = t4+times[5]
	if MPI.Initialized == true 
		times_add[6] = t5
	else
		times_add[6] = 0
	end
	times = times_add
    # for MPi version only
    #times_add[6] = t5
   ierr = 0 
   return niters, normr, normr0, ierr
end
