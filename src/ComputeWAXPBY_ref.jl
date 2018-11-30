#=
 @file ComputeWAXPBY_ref.cpp

 HPCG routine
=#

include("ComputeWAXPBY_ref.jl")
#=
  Routine to compute the update of a vector with the sum of two
  scaled vectors where: w = alpha*x + beta*y

  This is the reference WAXPBY impmentation.  It CANNOT be modified for the
  purposes of this benchmark.

  @param[in] n the number of vector elements (on this processor)
  @param[in] alpha, beta the scalars applied to x and y respectively.
  @param[in] x, y the input vectors
  @param[out] w the output vector.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeWAXPBY
=#
function ComputeWAXPBY_ref(n,alpha, x, beta, y, w) 

  @assert(length(x)>=n) # Test vector length
  @assert(length(y)>=n)

  xv = x
  yv = y
  wv = w

  if alpha==1.0
    for i=1:n
	wv[i] = xv[i] + beta * yv[i]
    end
  elif beta==1.0
    for i=1:n 
	wv[i] = alpha * xv[i] + yv[i]
    end
  else  
    for i=1:n 
	wv[i] = alpha * xv[i] + beta * yv[i]
    end
  end

  return 0
end
