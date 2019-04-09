#=
 @file ComputeWAXPBY_ref.cpp

 HPCG routine
=#

#=
  Routine to compute the update of a vector with the sum of two
  scaled vectors where: w = alpha*x + beta*y

  This is the reference WAXPBY impmentation.  It CANNOT be modified for the
  purposes of this benchmark.

  @param[out] w the output vector.
  @param[in] n the number of vector elements (on this processor)
  @param[in] alpha, beta the scalars applied to x and y respectively.
  @param[in] x, y the input vectors

  @return returns false upon success and true otherwise 

  @see compute_waxpby
=#

function compute_waxpby_ref!(w, n, alpha, x, beta, y)::Bool

  @assert(length(x)>=n) # Test vector length 
  @assert(length(y)>=n)

  if alpha == 1.0
    for i = 1:n
    	w[i] = x[i] + beta * y[i]
    end
  elseif beta == 1.0
    for i = 1:n 
	    w[i] = alpha * x[i] + y[i]
    end
  else  
    for i = 1:n 
	    w[i] = alpha * x[i] + beta * y[i]
    end
  end

  return false

end
