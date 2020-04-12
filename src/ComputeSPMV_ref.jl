
#=
  Precondition: First call exchange_externals to get off-processor values of x

  This is the reference SPMV implementation.  It CANNOT be modified for the
  purposes of this benchmark.

  @param[in]  A the known system matrix
  @param[in]  x the known vector
  @param[out] y the On exit contains the result: Ax.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeSPMV
  =#
 

  function compute_spmv_ref!(y::Array{Float64,1} , A::HPCGSparseMatrix, x::Array{Float64,1} ) 

  @assert(length(x)>=A.localNumberOfColumns) # Test vector lengths
  @assert(length(y)>=A.localNumberOfRows)

  if MPI.Initialized()== true
      exchange_halo!(x,A)
  end
  
#	  mv = Array{Float64,1}
#	  mIl  =Array{Int64,1}
  
  for i::Int64 = 1:A.localNumberOfRows
	  @inbounds mv = view(A.matrixValues, :, i)
	  @inbounds mIl = view(A.mtxIndL, :, i)
#	  nnz = Int64
@inbounds	  nnz   = A.nonzerosInRow[i]
	 for j::Int64= 1:nnz
		 @fastmath @inbounds y[i] +=  mv[j]*x[mIl[j]]
      end
  end
  return 0

end
