#= @file ComputeSYMGS_ref.cpp

 HPCG routine
=#

include("ExchangeHalo.jl")
include("ComputeSYMGS_ref.jl")

#=
Computes one step of symmetric Gauss-Seidel:

  Assumption about the structure of matrix A:
  - Each row 'i' of the matrix has nonzero diagonal value whose address is matrixDiagonal[i]
  - Entries in row 'i' are ordered such that:
       - lower triangular terms are stored before the diagonal element.
       - upper triangular terms are stored after the diagonal element.
       - No other assumptions are made about entry ordering.

  Symmetric Gauss-Seidel notes:
  - We use the input vector x as the RHS and start with an initial guess for y of all zeros.
  - We perform one forward sweep.  x should be initially zero on the first GS sweep, but we do not attempt to exploit this fact.
  - We then perform one back sweep.
  - For simplicity we include the diagonal contribution in the for-j loop, then correct the sum after

  @param[in] A the known system matrix
  @param[in] r the input vector
  @param[inout] x On entry, x should contain relevant values, on exit x contains the result of one symmetric GS sweep with r as the RHS.


  @warning Early versions of this kernel (Version 1.1 and earlier) had the r and x arguments in reverse order, and out of sync with other kernels.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeSYMGS
=#
function ComputeSYMGS_ref(A, r, x) 

  @assert(x.localLength==A.localNumberOfColumns) # Make sure x contain space for halo values

  ExchangeHalo(A,x)

  nrow = A.localNumberOfRows
  matrixDiagonal = A.matrixDiagonal  # An array of pointers to the diagonal entries A.matrixValues
  rv = r
  xv = x

  for i=1:nrow
    currentValues = A.matrixValues[i]
    currentColIndices = A.mtxIndL[i]
    currentNumberOfNonzeros = A.nonzerosInRow[i]
    currentDiagonal = matrixDiagonal[i][0] # Current diagonal value
    sum = rv[i] # RHS value

    for (int j=0 j< currentNumberOfNonzeros j++) 
      curCol = currentColIndices[j]
      sum -= currentValues[j] * xv[curCol]
    end
    sum += xv[i]*currentDiagonal # Remove diagonal contribution from previous loop

    xv[i] = sum/currentDiagonal

  end

  # Now the back sweep.

  for i=nrow-1: 0
    currentValues = A.matrixValues[i]
    currentColIndices = A.mtxIndL[i]
    currentNumberOfNonzeros = A.nonzerosInRow[i]
    currentDiagonal = matrixDiagonal[i][0] # Current diagonal value
    sum = rv[i] # RHS value

    for j = 1:currentNumberOfNonzeros
      curCol = currentColIndices[j]
      sum -= currentValues[j]*xv[curCol]
    end
    sum += xv[i]*currentDiagonal # Remove diagonal contribution from previous loop

    xv[i] = sum/currentDiagonal
  end

  return 0
end

