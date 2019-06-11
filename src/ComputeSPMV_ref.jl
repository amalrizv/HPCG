
using MPI
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
function compute_spmv_ref!(y, A, x) # takes SpMatrix_anx structure

  @assert(length(x)>=A.localNumberOfColumns) # Test vector lengths
  @assert(length(y)>=A.localNumberOfRows)

  if MPI.Initialized()
      exchange_halo!(x,A)
  end

  nrow = A.localNumberOfRows
  

  for i = 1:nrow
      sum      = 0
      cur_vals = A.matrixValues[i, :]
      cur_inds = A.mtxIndL[i, :]

      cur_nnz  = A.nonzerosInRow[i]

      for j= 1:cur_nnz


          sum = sum + (cur_vals[j]*x[cur_inds[j]])
      end
      y[i] = sum
  end
  
  return 0, y

end
