
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
function ComputeSPMV_ref(A, x, y) 

  @assert(length(x)>=A.localNumberOfColumns) # Test vector lengths
  @assert(length(y)>=A.localNumberOfRows)

  ExchangeHalo(A,x)
  xv = x.values
  yv = y.values
  nrow = A.localNumberOfRows
  for i=1:nrow
    sum = 0.0
    cur_vals = A.matrixValues[i]
    cur_inds = A.mtxIndL[i]
    cur_nnz = A.nonzerosInRow[i]
  end

    for j=1:cur_nnz
      sum += cur_vals[j]*xv[cur_inds[j]]
    end
    yv[i] = sum
  
  return 0
end
