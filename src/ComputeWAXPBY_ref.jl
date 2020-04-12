#=
 @file ComputeWAXPBY_ref.jl

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

function compute_waxpby_ref!(w::Array{Float64,1} , n::Int64, alpha::Float64, x::Array{Float64,1} , beta, y::Array{Float64,1} )

  @assert(length(x)>=n) # Test vector length 
  @assert(length(y)>=n)

  for i = 1:n 
  	if alpha == 1.0
	#	@show "alpha = 1.0"	  
	#	@show x[i], y[i]
    @fastmath @inbounds	w[i] =  x[i] + beta * y[i]

  	elseif beta == 1.0

	@fastmath @inbounds    w[i] = alpha * x[i] + y[i]
  	else  

	@fastmath @inbounds    w[i] = alpha * x[i] + beta * y[i]
    end

  end

  return 0

end
