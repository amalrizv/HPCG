
# Changelog
#
# Version 0.4
# - Added timing of setup time for sparse MV
# - Corrected percentages reported for sparse MV with overhead
#
####################################/

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
function test_cg!(A::HPCGSparseMatrix, data::CGData, b::Array{Float64, 1}, x::Array{Float64, 1},testcg_data::TestCGData ) 
    #testcg_data = TestCGData 
    # Use this array for collecting timing information
    times =  zeros(8)

    # Temporary storage for holding original diagonal and RHS
    exaggeratedDiagA = Array{Float64,1}(undef, A.localNumberOfRows)
    origB            = Array{Float64,1}(undef,A.localNumberOfRows)
    origDiagA        = Array{Float64,1}(undef, A.localNumberOfRows)
    fill!(origDiagA, 26.0)
    exaggeratedDiagA[1:length(origDiagA)] = origDiagA
    origB[1:length(b)]           = b
       # Modify the matrix diagonal to greatly exaggerate diagonal values.
    # CG should converge in about 10 iterations for this problem, regardless of problem size
    for i=1:A.localNumberOfRows
        globalRowID = A.localToGlobalMap[i]
        if globalRowID<10
            scale = (globalRowID+1)*1.0e6
			exaggeratedDiagA[i] = exaggeratedDiagA[i] .* scale
			b[i] = b[i] .* scale
        else 
			exaggeratedDiagA[i]  = exaggeratedDiagA[i] .* 1.0e6
			b[i] = b[i] .* 1.0e6
        end
    end

    for i = 1: A.localNumberOfRows 
	A.matrixValues[A.curcols[i], i] = exaggeratedDiagA[i]
    end

    niters          = 0
    normr           = 0.0
    normr0          = 0.0
    maxIters        = 50
    numberOfCgCalls = 2
    tolerance       = 1.0e-12 # Set tolerance to reasonable value for grossly scaled diagonal terms
    
    # For the unpreconditioned CG call, we should take about 10 iterations, permit 12
    # For the preconditioned case, we should take about 1 iteration, permit 2
    
    for k=1:2 # This loop tests both unpreconditioned and preconditioned runs
        expected_niters = testcg_data.expected_niters_no_prec
        if k==2
            expected_niters = testcg_data.expected_niters_prec
        end
        for i=1:numberOfCgCalls
			zero_fill!(x)
			# CG is being called 4 times 
        	niters, normr, normr0, ierr = cg!(A, data, b, x, maxIters, tolerance, times, k==2)
            if ierr==1
                @debug("Error in call to CG:$ierr.\n")
            end
            if niters <= expected_niters 
                testcg_data.count_pass+=1
            else 
                testcg_data.count_fail+=1
            end
            if k==1 && niters>testcg_data.niters_max_no_prec
                testcg_data.niters_max_no_prec = niters # Keep track of largest iter count
            end
            if k==2 && niters>testcg_data.niters_max_prec
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
    for i = 1: A.localNumberOfRows 
	A.matrixValues[A.curcols[i], i] = exaggeratedDiagA[i]
    end

	b[1:length(origB)] = origB

    # Delete vectors
    origDiagA         = nothing
    exaggeratedDiagA  = nothing
    origB             = nothing
    testcg_data.normr = normr
	return testcg_data

end
