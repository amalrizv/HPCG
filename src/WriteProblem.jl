#=
  Routine to dump:
   - matrix in row, col, val format for analysis with MATLAB
   - x, xexact, b as simple arrays of numbers.

   Writes to A.dat, x.dat, xexact.dat and b.dat, respectivly.

   NOTE:  THIS CODE ONLY WORKS ON SINGLE PROCESSOR RUNS

   Read into MATLAB using:

       load A.dat
       A=spconvert(A)
       load x.dat
       load xexact.dat
       load b.dat

  @param[in] geom   The description of the problem's geometry.
  @param[in] A      The known system matrix
  @param[in] b      The known right hand side vector
  @param[in] x      The solution vector computed by CG iteration
  @param[in] xexact Generated exact solution

  @return Returns with -1 if used with more than one MPI process. Returns with 0 otherwise.

  @see GenerateProblem
=#

function WriteProblem(geom, A, b, x, xexact) 

  if geom.size!=1 
	return -1 #TODO Only works on one processor.  Need better error handler
  end
  nrow = A.totalNumberOfRows

  fx = fopen("x.dat", "w")
  fxexact = fopen("xexact.dat", "w")
  fb = fopen("b.dat", "w")

  if ! fA || ! fx || ! fxexact || ! fb 
    if fb 
	fclose(fb)
    end
    if fxexact
	fclose(fxexact)
    end
    if fx
	fclose(fx)
    end
    if fA 
	fclose(fA)
    end
    return -1
  end

  for i=1:nrow
    currentRowValues = A.matrixValues[i,:]
    currentRowIndices = A.mtxIndG[i,:]
    currentNumberOfNonzeros = A.nonzerosInRow[i]
    for j=1:currentNumberOfNonzeros
#ifdef HPCG_NO_LONG_LONG
      fprintf(fA, " %d %d %22.16e\n",i+1,(global_int_t)(currentRowIndices[j]+1),currentRowValues[j])
#else
      fprintf(fA, " %lld %lld %22.16e\n",i+1,(global_int_t)(currentRowIndices[j]+1),currentRowValues[j])
#endif
    end
    fprintf(fx, "%22.16e\n",x.values[i])
    fprintf(fxexact, "%22.16e\n",xexact.values[i])
    fprintf(fb, "%22.16e\n",b.values[i])
    end

  fclose(fA)
  fclose(fx)
  fclose(fxexact)
  fclose(fb)
  return 0
end
