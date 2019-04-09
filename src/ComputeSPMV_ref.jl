
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
function ComputeSPMV_ref(AAA, x, y) #takes SpMatrix_anx structure

  AA = AAA.sp_matrix
  A = AA.sp_matrix
  #@assert(length(x)>=AAA.localNumberOfCols) # Test vector lengths
  #@assert(length(y)>=AA.localNumberOfRows)

  ExchangeHalo(AAA,x)
  xv = x
  yv = y
  nrow = AA.localNumberOfRows
  for i=1:nrow
    sum = 0
    cur_vals = AA.matrixValues[i]
    cur_inds = AA.mtxIndL[i]
    cur_nnz = AA.nonzerosInRow[i]

    for j=1:cur_nnz
       
#      sum = sum+(cur_vals[j]*xv[cur_inds[j]])
    end
#    yv[i] = sum
  end
  
  return 0
end
