
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
function compute_spmv_ref!(y, A, x) 
  @assert(length(x)>=A.localNumberOfColumns) # Test vector lengths
  @assert(length(y)>=A.localNumberOfRows)

  if A.geom.rank == 0
	  fs =  open("mpi_ircv_0.txt", "a")
	  fp =  open("mpi_send_0.txt", "a")
  else
 	  fs =  open("mpi_ircv_1.txt", "a")
	  fp =  open("mpi_send_1.txt", "a")
  end
  if MPI.Initialized()== true
	  println(fs, "Called from : SPMV")
	  println(fp, "Called from : SPMV")
      exchange_halo!(x,A)
  end

  nrow = A.localNumberOfRows
  

  for i = 1:nrow
      sum      = 0
      cur_vals = A.matrixValues[i, :]
      cur_inds = A.mtxIndL[i, :]

      cur_nnz  = A.nonzerosInRow[i]
	  # print output of both ranks 


      for j= 1:cur_nnz
#	  if A.geom.rank == 0
#		  open("spmv_0.txt", "a") do f
#			  println(f, "sum = $sum + cur_vals[$j]( = $(cur_vals[j]) ) * x[cur_inds[$j] ( = x[$(cur_inds[j])] = $(x[cur_inds[j]]) )")
#		  end
#	  else
#		  open("spmv_1.txt", "a") do f
#			  println(f, "sum = $sum + cur_vals[$j]( = $(cur_vals[j]) ) * x[cur_inds[$j] ( = x[$(cur_inds[j])] = $(x[cur_inds[j]]) )")
#		  end
#	  end

          sum = sum + (cur_vals[j]*x[cur_inds[j]])

      end
      y[i] = sum
  end
#		if A.geom.rank == 0 
#			open("b_computed_0.txt", "a") do f
#				println(f,"spmv_b_computed[$(length(y))] = $(y[length(y)])")
#			end
#		else
#			open("b_computed_1.txt", "a") do f
#				println(f,"spmv_b_computed[$(length(y))] = $(y[length(y)])")
#			end
#		end

  return 0

end
