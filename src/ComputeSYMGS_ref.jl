#= @file ComputeSYMGS_ref.cpp

 HPCG routine
=#

include("ExchangeHalo.jl")


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
  @param[inout] xv On entry, xv should contain relevant values, on exit xv contains the result of one symmetric GS sweep with r as the RHS.


  @warning Early versions of this kernel (Version 1.1 and earlier) had the r and x arguments in reverse order, and out of sync with other kernels.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeSYMGS
=#
function compute_symgs_ref!(xv, A, rv) 
 # xv = x
  @assert(length(xv)==A.localNumberOfColumns) # Make sure x contain space for halo values
  if A.geom.rank == 0
	  fs =  open("mpi_ircv_0.txt", "a")
	  fp =  open("mpi_send_0.txt", "a")
	  sfs =  open("symgs_front_swp_0.txt", "a")
  else
 	  fs =  open("mpi_ircv_1.txt", "a")
	  fp =  open("mpi_send_1.txt", "a")
	  sfs =  open("symgs_front_swp_1.txt", "a")
  end
  if MPI.Initialized()== true
	  println(fs, "Called from SYMGS")
	  println(fp, "Called from SYMGS")
	  exchange_halo!(xv,A)
  end
  close(fs)
  close(fp)

  nrow = A.localNumberOfRows
  #matrixDiagonal = A.matrixDiagonal  # An array of pointers to the diagonal entries A.matrixValues
  matrixValues = A.matrixValues
  mtxIndL = A.mtxIndL
	  println(sfs, "xv_original[1] = $(xv[1])")
	  println(sfs, "sum = sum - currentValues[j] * xv[curCol]")

  for i=1:nrow
    currentValues = matrixValues[i, :]
    currentColIndices = mtxIndL[i, :]
    currentNumberOfNonzeros = A.nonzerosInRow[i]
    curcols = A.curcols
    currentDiagonal = matrixValues[i,curcols[i] ] # Current diagonal value
    sum = rv[i] # RHS value
#	println(sfs, "sum = sum - currentValues[j] * xv[curCol]")
    for j=1:currentNumberOfNonzeros 
      curCol = currentColIndices[j]
	  println(sfs, "sum = $sum - $(currentValues[j]) * $(xv[curCol])")
	  # RZV First iteration of this loop has a different value 
      sum = sum - currentValues[j] * xv[curCol]
    end
#	println(sfs, " sum =sum + xv[i]*currentDiagonal")
#	println(sfs, "sum = $sum + $(xv[i]) * $currentDiagonal")
    sum =sum + xv[i]*currentDiagonal # Remove diagonal contribution from previous loop
#	println(sfs, "xv[i] = sum/currentDiagonal")
#	println(sfs, "xv[i] = $sum / $currentDiagonal")
    xv[i] = sum/currentDiagonal

  end
#RZV

#	  println(sfs, "xv_modified[1] = $(xv[1])")
		
  # Now the back sweep.

  for i=Iterators.reverse(1:nrow)
    currentValues = matrixValues[i, :]
    currentColIndices = mtxIndL[i,:]
    currentNumberOfNonzeros = A.nonzerosInRow[i]
    curcols = A.curcols
    currentDiagonal = matrixValues[i,curcols[i] ] # Current diagonal value
    sum = rv[i] # RHS value
#	println(sfs, "sum = sum - currentValues[j] * xv[curCol]")

    for j = 1:currentNumberOfNonzeros
      curCol = currentColIndices[j]
#	  println(sfs, "sum = $sum - $(currentValues[j]) * $(xv[curCol])")
      sum = sum - currentValues[j]*xv[curCol]
    end
#	println(sfs, " sum =sum + xv[i]*currentDiagonal")
#    println(sfs, "sum = $sum + $(xv[i]) * $currentDiagonal")
    sum =sum + xv[i]*currentDiagonal # Remove diagonal contribution from previous loop
#	println(sfs, " sum =sum + xv[i]*currentDiagonal")
#	println(sfs, "xv[i] = $sum / $currentDiagonal")
    xv[i] = sum/currentDiagonal

  end
#	  println(sfs, "xv_back_modified[1] = $(xv[1])")
	close(sfs)
  return 0
end

