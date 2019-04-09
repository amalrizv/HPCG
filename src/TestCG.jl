
# Changelog
#
# Version 0.4
# - Added timing of setup time for sparse MV
# - Corrected percentages reported for sparse MV with overhead
#
####################################/

include("hpcg.jl")

include("TestCG_struct.jl")
include("CG.jl")

#=
  Test the correctness of the Preconditined CG implementation by using a system matrix with a dominant diagonal.

  @param[in]    geom The description of the problem's geometry.
  @param[in]    A    The known system matrix
  @param[in]    data the data structure with all necessary CG vectors preallocated
  @param[in]    b    The known right hand side vector
  @param[inout] x    On entry: the initial guess on exit: the new approximate solution
  @param[out]   testcg_data the data structure with the results of the test including pass/fail information

  @return Returns zero on success and a non-zero value otherwise.

  @see CG()
=#
function TestCG(AAAA, data, b, x, count_pass, count_fail) 
  println("here")
  println(typeof(AAAA))
  AAA = AAAA.sp_matrix #sp_anx structure
  AA = AAA.sp_matrix   #sp_matrix structure
  A = AA.sp_matrix     #sp_init structure
  testcg_data=TestCGData 
  # Use this array for collecting timing information
  times =  zeros(8)
  # Temporary storage for holding original diagonal and RHS
  exaggeratedDiagA = Vector{Int64}(undef, AA.localNumberOfRows)
  origB =Vector(undef,AA.localNumberOfRows)
  origDiagA = CopyMatrixDiagonal(AA)
  exaggeratedDiagA = origDiagA
  origB = b

  # Modify the matrix diagonal to greatly exaggerate diagonal values.
  # CG should converge in about 10 iterations for this problem, regardless of problem size
  for i=1:AA.localNumberOfRows
    globalRowID = AA.localToGlobalMap[i]
    if globalRowID<9
      scale = (globalRowID+2)*1.0e6
      exaggeratedDiagA = exaggeratedDiagA .* scale
      b = b .* scale
     else 
      exaggeratedDiagA  = exaggeratedDiagA .* 1.0e6
      b = b .* 1.0e6
    end
  end
  AA = ReplaceMatrixDiagonal(AA, exaggeratedDiagA)

  niters = 0
  normr = 0.0
  normr0 = 0.0
  maxIters = 50
  numberOfCgCalls = 2
  tolerance = 1.0e-12 # Set tolerance to reasonable value for grossly scaled diagonal terms
  testcg_data = TestCGData(count_pass, count_fail,12,2,0,0, normr) 
  
  # For the unpreconditioned CG call, we should take about 10 iterations, permit 12
  # For the preconditioned case, we should take about 1 iteration, permit 2
  
  for k=0:2 # This loop tests both unpreconditioned and preconditioned runs
    expected_niters = testcg_data.expected_niters_no_prec
    if k==1
	 expected_niters = testcg_data.expected_niters_prec
    end
    for i=0:numberOfCgCalls
      x = zeros(length(x)) # Zero out x
      ierr, times_add = CG(AAAA, data, b, x, maxIters, tolerance, niters, normr, normr0, times, k==1)
      times = times_add
      if ierr==1
	@debug("Error in call to CG:$ierr.\n")
      end
      if niters <= expected_niters 
        testcg_data.count_pass+=1
       else 
        testcg_data.count_fail+=1
      end
      if k==0 && niters>testcg_data.niters_max_no_prec
	 testcg_data.niters_max_no_prec = niters # Keep track of largest iter count
      end
      if k==1 && niters>testcg_data.niters_max_prec
	 testcg_data.niters_max_prec = niters # Same for preconditioned run
      end
      if A.geom.rank==0
        @debug("Call [$i] Number of Iterations [$niters] Scaled Residual [$(normr/normr0)]")
        if niters > expected_niters
          @debug(" Expected $expected_niters iterations.  Performed $niters .")
        end
      end
    end
  end

  # Restore matrix diagonal and RHS
  ReplaceMatrixDiagonal(A, origDiagA)
  b = origB
  # Delete vectors
  origDiagA = nothing
  exaggeratedDiagA= nothing
  origB = nothing
  testcg_data.normr = normr

  return testcg_data, times
end
