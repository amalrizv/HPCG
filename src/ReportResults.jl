#=
 @file ReportResults.cpp

 HPCG routine
=#

using Dates

#=
 Creates a YAML file and writes the information about the HPCG run, its results, and validity.

  @param[in] geom The description of the problem's geometry.
  @param[in] A    The known system matrix
  @param[in] numberOfMgLevels Number of levels in multigrid V cycle
  @param[in] numberOfCgSets Number of CG runs performed
  @param[in] niters Number of preconditioned CG iterations performed to lower the residual below a threshold
  @param[in] times  Vector of cumulative timings for each of the phases of a preconditioned CG iteration
  @param[in] testcg_data    the data structure with the results of the CG-correctness test including pass/fail information
  @param[in] testsymmetry_data the data structure with the results of the CG symmetry test including pass/fail information
  @param[in] testnorms_data the data structure with the results of the CG norm test including pass/fail information
  @param[in] global_failure indicates whether a failure occurred during the correctness tests of CG

  @see YAML_Doc
=#
function report_results(A, numberOfMgLevels, numberOfCgSets, refMaxIters, optMaxIters, times, testcg_data, testsymmetry_data, testnorms_data, global_failure, quickPath) 
	  

  minOfficialTime = 1800 # Any official benchmark result must run at least this many seconds

  if MPI.Initialized() == true
 		t4 = times[4]
 		t4min = 0.0
 		t4max = 0.0
		t4avg = 0.0
  		t4min 	= MPI.Allreduce(t4,  MPI.MIN, MPI.COMM_WORLD)
  		t4max 	= MPI.Allreduce(t4,  MPI.MAX, MPI.COMM_WORLD)
  		t4avg		= MPI.Allreduce(t4,  MPI.SUM, MPI.COMM_WORLD)
  		t4avg = t4avg/A.geom.size
 	end

 	if A.geom.rank==0  # Only PE 0 needs to compute and report timing results

    # TODO: Put the FLOP count, Memory BW and Memory Usage models into separate functions
	#		 ======================== FLOP count model =======================================

    	fNumberOfCgSets 	= numberOfCgSets
    	fniters 			= fNumberOfCgSets *  optMaxIters
    	fnrow 			= A.totalNumberOfRows
    	fnnz 			= A.totalNumberOfNonzeros

   # Op counts come from implementation of CG in CG.cpp (include 1 extra for the CG preamble ops)

    	fnops_ddot 		= (3.0*fniters+fNumberOfCgSets)*2.0*fnrow # 3 ddots with nrow adds and nrow mults
    	fnops_waxpby 	= (3.0*fniters+fNumberOfCgSets)*2.0*fnrow # 3 WAXPBYs with nrow adds and nrow mults
    	fnops_sparsemv 	= (fniters+fNumberOfCgSets)*2.0*fnnz # # SpMV with nnz adds and nnz mults

   # Op counts from the multigrid preconditioners

    	fnops_precond = 0.0
    	Af = A

    	for i=1:numberOfMgLevels-1
      		fnnz_Af 					= Af.totalNumberOfNonzeros
      		fnumberOfPresmootherSteps 	= Af.mgData.numberOfPresmootherSteps
      		fnumberOfPostsmootherSteps = Af.mgData.numberOfPostsmootherSteps
      		fnops_precond += fnumberOfPresmootherSteps*fniters*4.0*fnnz_Af # number of presmoother flops
      		fnops_precond += fniters*2.0*fnnz_Af # cost of fine grid residual calculation
      		fnops_precond += fnumberOfPostsmootherSteps*fniters*4.0*fnnz_Af  # number of postsmoother flops
      		Af = Af.Ac 	# Go to next coarse level
    	end

    	fnops_precond 	+= fniters*4.0*( Af.totalNumberOfNonzeros) # One symmetric GS sweep at the coarsest level
    	fnops 			 = fnops_ddot+fnops_waxpby+fnops_sparsemv+fnops_precond
    	frefnops 		 = fnops * ( refMaxIters)/( optMaxIters)

    # ======================== Memory bandwidth model =======================================

    # Read/Write counts come from implementation of CG in CG.cpp (include 1 extra for the CG preamble ops)
	
    	fnreads_ddot 	= (3.0*fniters+fNumberOfCgSets)*2.0*fnrow*sizeof(Float64) # 3 ddots with 2 nrow reads
    	fnwrites_ddot 	= (3.0*fniters+fNumberOfCgSets)*sizeof(Float64) # 3 ddots with 1 write
    	fnreads_waxpby 	= (3.0*fniters+fNumberOfCgSets)*2.0*fnrow*sizeof(Float64) # 3 WAXPBYs with nrow adds and nrow mults
    	fnwrites_waxpby 	= (3.0*fniters+fNumberOfCgSets)*fnrow*sizeof(Float64)# 3 WAXPBYs with nrow adds and nrow mults
    	fnreads_sparsemv = (fniters+fNumberOfCgSets)*(fnnz*(sizeof(Float64)+sizeof(Int64)) + fnrow*sizeof(Float64))# 1 SpMV with nnz reads of values, nnz reads indices,

    	# plus nrow reads of x
    	fnwrites_sparsemv = (fniters+fNumberOfCgSets)*fnrow*sizeof(Float64) # 1 SpMV nrow writes
    	# Op counts from the multigrid preconditioners
    	fnreads_precond 	= 0.0
    	fnwrites_precond	= 0.0
   		Af = A
    	for i=1:numberOfMgLevels-1 
     		fnnz_Af 		= Af.totalNumberOfNonzeros
     		fnrow_Af 		= Af.totalNumberOfRows
     		fnumberOfPresmootherSteps 	= Af.mgData.numberOfPresmootherSteps
     		fnumberOfPostsmootherSteps = Af.mgData.numberOfPostsmootherSteps
     		fnreads_precond 	+= fnumberOfPresmootherSteps*fniters*(2.0*fnnz_Af*(sizeof(Float64)+sizeof(Int64)) + fnrow_Af*sizeof(Float64)) # number of presmoother reads
     		fnwrites_precond 	+= fnumberOfPresmootherSteps*fniters*fnrow_Af*sizeof(Float64) # number of presmoother writes
     		fnreads_precond 	+= fniters*(fnnz_Af*(sizeof(Float64)+sizeof(Int64)) + fnrow_Af*sizeof(Float64)) # Number of reads for fine grid residual calculation
     		fnwrites_precond 	+= fniters*fnnz_Af*sizeof(Float64) # Number of writes for fine grid residual calculation
     		fnreads_precond 	+= fnumberOfPostsmootherSteps*fniters*(2.0*fnnz_Af*(sizeof(Float64)+sizeof(Int64)) + fnrow_Af*sizeof(Float64))  # number of postsmoother reads
     		fnwrites_precond 	+= fnumberOfPostsmootherSteps*fniters*fnnz_Af*sizeof(Float64)  # number of postsmoother writes
      		Af = Af.Ac # Go to next coarse level
  		end

    	fnnz_Af 	= Af.totalNumberOfNonzeros
    	fnrow_Af 	= Af.totalNumberOfRows
    	fnreads_precond += fniters*(2.0*fnnz_Af*(sizeof(Float64)+sizeof(Int64)) + fnrow_Af*sizeof(Float64)) # One symmetric GS sweep at the coarsest level
    	fnwrites_precond += fniters*fnrow_Af*sizeof(Float64) # One symmetric GS sweep at the coarsest level
    	fnreads 	= fnreads_ddot+fnreads_waxpby+fnreads_sparsemv+fnreads_precond
    	fnwrites 	= fnwrites_ddot+fnwrites_waxpby+fnwrites_sparsemv+fnwrites_precond
    	frefnreads 	= fnreads * ( refMaxIters)/( optMaxIters)
    	frefnwrites = fnwrites * ( refMaxIters)/( optMaxIters)


    # ======================== Memory usage model =======================================

    # Data in GenerateProblem_ref

   		numberOfNonzerosPerRow = 27.0 # We are approximating a 27-point finite element/volume/difference 3D stencil
    	size = A.geom.size # Needed for estimating size of halo

    	fnbytes 		= ( sizeof(Geometry))      # Geometry struct in main.cpp
    	fnbytes 	+= ( sizeof(Float64)*fNumberOfCgSets) # testnorms_data in main.cpp

    # Model for GenerateProblem_ref.cpp
    	fnbytes += fnrow*sizeof(Int64)      # array nonzerosInRow
    	fnbytes += fnrow*( sizeof(Int64)) # mtxIndG
    	fnbytes += fnrow*( sizeof(Int64))  # mtxIndL
    	fnbytes += fnrow*( sizeof(Float64))      # matrixValues
    	fnbytes += fnrow*( sizeof(Float64))      # matrixDiagonal
    	fnbytes += fnrow*numberOfNonzerosPerRow*( sizeof(Int64))  # mtxIndL[1..nrows]
    	fnbytes += fnrow*numberOfNonzerosPerRow*( sizeof(Float64))       # matrixValues[1..nrows]
    	fnbytes += fnrow*numberOfNonzerosPerRow*( sizeof(Int64)) # mtxIndG[1..nrows]
    	fnbytes += fnrow*( 3*sizeof(Float64)) # x, b, xexact

    # Model for CGData.hpp
   		fncol = ( A.localNumberOfColumns) * size # Estimate of the global number of columns using the value from rank 0
    	fnbytes += fnrow*( 2*sizeof(Float64)) # r, Ap
    	fnbytes += fncol*( 2*sizeof(Float64)) # z, p

		fnbytesPerLevel  	= Vector{Float64}(undef, numberOfMgLevels) # Count byte usage per level (level 0 is main CG level)
    	fnbytesPerLevel[1] 	= fnbytes

    # Benchmarker-provided model for OptimizeProblem.cpp
    	fnbytes_OptimizedProblem = OptimizeProblemMemoryUse(A)
    	fnbytes += fnbytes_OptimizedProblem

    	Af = A.Ac
    	for i=1:numberOfMgLevels-1 
     		fnrow_Af = Af.totalNumberOfRows
     		fncol_Af = ( Af.localNumberOfColumns) * size # Estimate of the global number of columns using the value from rank 0
     		fnbytes_Af = 0.0
      # Model for GenerateCoarseProblem.cpp
     		fnbytes_Af += fnrow_Af*( sizeof(Int64)) # f2cOperator
    		fnbytes_Af += fnrow_Af*( sizeof(Float64)) # rc
     		fnbytes_Af += 2.0*fncol_Af*( sizeof(Float64)) # xc, Axf are estimated based on the size of these arrays on rank 0
     		fnbytes_Af += sizeof(Geometry)+sizeof(HPCGSparseMatrix)#=+3*sizeof(Vector)=#+sizeof(MGData) # Account for structs geomc, Ac, rc, xc, Axf - (minor)

      # Model for GenerateProblem.cpp (called within GenerateCoarseProblem.cpp)
      		fnbytes_Af += fnrow_Af*sizeof(Int64)      # array nonzerosInRow
      		fnbytes_Af += fnrow_Af*( sizeof(Int64)) # mtxIndG
      		fnbytes_Af += fnrow_Af*( sizeof(Int64))  # mtxIndL
      		fnbytes_Af += fnrow_Af*( sizeof(Float64))      # matrixValues
      		fnbytes_Af += fnrow_Af*( sizeof(Float64))      # matrixDiagonal
      		fnbytes_Af += fnrow_Af*numberOfNonzerosPerRow*( sizeof(Int64))  # mtxIndL[1..nrows]
      		fnbytes_Af += fnrow_Af*numberOfNonzerosPerRow*( sizeof(Float64))       # matrixValues[1..nrows]
      		fnbytes_Af += fnrow_Af*numberOfNonzerosPerRow*( sizeof(Int64)) # mtxIndG[1..nrows]

      # Model for SetupHalo_ref.cpp
	  		if MPI.Initialized() ==  true 
       			fnbytes_Af += ( sizeof(Float64)*Af.totalToBeSent) #sendBuffer
       			fnbytes_Af += ( sizeof(Int64)*Af.totalToBeSent) # elementsToSend
       			fnbytes_Af += ( sizeof(Int64)*Af.numberOfSendNeighbors) # neighbors
       			fnbytes_Af += ( sizeof(Int64)*Af.numberOfSendNeighbors) # receiveLength, sendLength
	  		end
      		fnbytesPerLevel[i] = fnbytes_Af
      		fnbytes	+= fnbytes_Af # Running sum
      		Af 		= Af.Ac # Go to next coarse level
    	end

    	#@assert(Af==0) # Make sure we got to the lowest grid level

    # Count number of bytes used per equation
    	fnbytesPerEquation = fnbytes/fnrow

    # Instantiate YAML document
		date=  now()
		report_result_filename  ="HPCG-Benchmark-Julia_"*string(date)*".txt"
	  	report_result = open(report_result_filename, "w")
		hpcg_specs  = "HPCG-Benchmark\n HPCG-Benchmark, 3.1\n Release date\n"
		

    	machine_summary = "Machine Summary\n\tDistributed Processes $(A.geom.size)\n\tThreads per processes,$(A.geom.numThreads)\n"

	    global_problem_dimensions= "Global Problem Dimensions \n \tGlobal nx,$(A.geom.gnx)\n\tGlobal ny, $(A.geom.gny)\n\tGlobal nz, $(A.geom.gnz)\n"

    	processor_dimensions = "Processor Dimensions \n\tnpx $(A.geom.npx)\n\tnpy, $(A.geom.npy)\n\tnpz, $(A.geom.npz)\n"

    	local_domain_dimensions = "Local Domain Dimensions\n\tnx ,$(A.geom.nx) \n\tny ,$(A.geom.ny)\n"
		println(report_result, hpcg_specs, machine_summary, global_problem_dimensions, processor_dimensions, local_domain_dimensions)
    	ipartz_ids = 1

    for i=1:A.geom.npartz
      println(report_result, "Lower ipz", ipartz_ids)
      println(report_result, "Upper ipz", A.geom.partz_ids[i]-1)
      println(report_result, "\tnz ",A.geom.partz_nz[i])
      ipartz_ids = A.geom.partz_ids[i]
    end

    	title_problem_summary = "########## Problem Summary  ##########\n"

    	setup_info = "Setup Information\n\tSetup Time $(times[10])\n"

    	linear_system_information = "Linear System Information\n\tNumber of Equations $(A.totalNumberOfRows) \n\tNumber of Nonzero Terms , $(A.totalNumberOfNonzeros)\n"
    	Af = A

    	mg_info  = "Multigrid Information\n\tNumber of coarse grid levels, $(numberOfMgLevels-1) \n Coarse Grids\n"
		println(report_result, title_problem_summary, setup_info, linear_system_information, mg_info)
	
    for i=1:numberOfMgLevels-1
      println(report_result, "\tGrid Level $i")
	  println(report_result, "\tNumber of Equations $(Af.Ac.totalNumberOfRows)")
	  println(report_result, "\tNumber of Nonzero Terms $(Af.Ac.totalNumberOfNonzeros)")
	  println(report_result, "\tNumber of Presmoother Steps $(Af.mgData.numberOfPresmootherSteps)")
	  println(report_result, "\tNumber of Postsmoother Steps $(Af.mgData.numberOfPostsmootherSteps)")
      Af = Af.Ac
    end
    	title_memory_use = "########## Memory Use Summary  ##########\n"

    	mem_use_info = "Memory Use Information\n\tTotal memory used for data (Gbytes),$(fnbytes/1000000000.0)\n\tMemory used for OptimizeProblem data (Gbytes), $(fnbytes_OptimizedProblem/1000000000.0)\n\tBytes per equation (Total memory / Number of Equations),$(fnbytesPerEquation) Memory used for linear system and CG (Gbytes), $(fnbytesPerLevel[1]/1000000000.0)\n\tCoarse Grids\n"
		println(report_result, title_memory_use, mem_use_info)

    for i=1:numberOfMgLevels-1
      println(report_result, "\tGrid Level $i")
	  println(report_result, "\tMemory used $(fnbytesPerLevel[i]/1000000000.0)")
    end

    	title_v_v_testing  = "########## V&V Testing Summary  ##########\n"

		if (testcg_data.count_fail==0)
			test_cg_result = "PASSED"
    	else
			test_cg_result = "FAILED"
		end

    	spectral_convergence = "Spectral Convergence Tests\t"* test_cg_result *"\n Unpreconditioned\n\tMaximum iteration count $( testcg_data.niters_max_no_prec)\n\tExpected iteration coun $( testcg_data.expected_niters_no_prec)\n Preconditioned \n\tMaximum iteration count $(testcg_data.niters_max_prec) \n\tExpected iteration count $( testcg_data.expected_niters_prec)\n"

#@show testsymmetry_data
		if (testsymmetry_data.count_fail==0)
			test_symmetry_result = "PASSED"
    	else
			test_symmetry_result = "FAILED"
		end

   
		departure_from_symmetry = "Departure from Symmetry |x'Ay-y'Ax|/(2*||x||*||A||*||y||)/epsilon\t"*test_symmetry_result*"\n\tDeparture for SpMV $( testsymmetry_data.depsym_spmv)\n\tDeparture for MG $(testsymmetry_data.depsym_mg)\n"

	    title_iteration_summary = "########## Iterations Summary  ##########\n"
 
		if global_failure == 0 
			iterations_result = "PASSED"
    	else
			iterations_result = "FAILED"
		end
     	iteration_count_information = "Iteration Count Information\t" *iterations_result* "\n\tReference CG iterations per set$( refMaxIters) \n \tOptimized CG iterations per set $(optMaxIters)\n \t Total number of reference iterations $(refMaxIters*numberOfCgSets)\n\tTotal number of optimized iterations $( optMaxIters*numberOfCgSets)\n"

    	title_reproducibility = "########## Reproducibility Summary  ##########\n"
	 	if testnorms_data.pass ==1
      		testnorms_result =  "PASSED"
    	else
     	 	testnorms_result =  "FAILED"
		end
		reproducibility_information = "Reproducibility Information\t"*testnorms_result *"\n\tScaled residual mean $(testnorms_data.mean) \n\tScaled residual variance $(testnorms_data.variance)\n"

    	title_performance_summary = "########## Performance Summary (times in sec) ##########\n"

    	benchmark_time_summary = " Benchmark Time Summary\n\tOptimization phase $(times[8]) \n\tDDOT $(times[2])\n\tWAXPBY $(times[3])\n\tSpMV $(times[4])\n\tMG $(times[6])\n\tTotal $(times[1])\n"

    	floating_point_ops_summary = " Floating Point Operations Summary \n\tRaw DDOT $(fnops_ddot)\n\tRaw WAXPBY$(fnops_waxpby)\n\tRaw SpMV $(fnops_sparsemv)\n\tRaw MG $(fnops_precond)\n\tTotal $(fnops)\n\tTotal with convergence overhead $(frefnops)\n"

	    total = (frefnreads+frefnwrites)/(times[1]+fNumberOfCgSets*(times[7]/10.0+times[9]/10.0))/1.0e9
		gb_s_summary = "GB/s Summary \n\tRaw Read B/W $(fnreads/times[1]/1.0E9)\n\tRaw Write B/W $(fnwrites/times[1]/1.0e9)\n\tRaw Total B/W $((fnreads+fnwrites)/(times[1])/1.0E9))\n\tTotal with convergence and optimization phase overhead $total \n"

 # This final GFLOP/s rating includes the overhead of
 # problem setup and optimizing the data structures 
 # vs ten sets of 50 iterations of CG
 
   		totalGflops = frefnops/(times[1]+fNumberOfCgSets*(times[8]/10.0+times[8]/10.0))/1.0E9
   		totalGflops24 = frefnops/(times[1]+fNumberOfCgSets*times[8]/10.0)/1.0E9


		g_flops_summary = "GFLOP/s Summary\n\tRaw DDOT $(fnops_ddot/times[2]/1.0E9)\n\tRaw WAXPBY $(fnops_waxpby/times[3]/1.0E9)\n\tRaw SpMV $(fnops_sparsemv/(times[4])/1.0E9)\n\tRaw MG $(fnops_precond/(times[6])/1.0E9)\n\tRaw Total $(fnops/times[1]/1.0E9)\n\tTotal with convergence overhead $(frefnops/times[1]/1.0E9)\n\tTotal with convergence and optimization phase overhead $(totalGflops)\n"

    	user_opt_info = "User Optimization Overheads\n\tOptimization phase time (sec)$(times[7])\n\tOptimization phase time vs reference SpMV+MG time $(times[7]/times[8])\n"
		 println(report_result, title_v_v_testing,  spectral_convergence,  departure_from_symmetry, title_iteration_summary, iteration_count_information, title_reproducibility,  reproducibility_information, title_performance_summary, benchmark_time_summary, floating_point_ops_summary ,gb_s_summary, g_flops_summary,  user_opt_info)
		if MPI.Initialized() == true
			ddot_timing_variations = "DDOT Timing Variations\n\tMin DDOT MPI_Allreduce time $(t4min)\n\tMax DDOT MPI_Allreduce time $(t4max) \n\tAvg DDOT MPI_Allreduce time $(t4avg)\n\t"

    		#sparse_operations_overheads = "Sparse Operations Overheads\n"
			#halo_exchange_time = "\tHalo exchange time (sec) $(times[7])\n"
    		#halo_exchange_spmv_percentage = "\tHalo exchange as percentage of SpMV time $((times[7])/totalSparseMVTime*100.0)\n"
			println(report_result, ddot_timing_variations)
			#println(report_result, sparse_operations_overheads , halo_exchange_time, halo_exchange_spmv_percentage)
		end
     	final_summary_title = "Final Summary\n"

    	isValidRun = (testcg_data.count_fail==0) && (testsymmetry_data.count_fail==0) && (testnorms_data.pass ==1) && (global_failure ==0)
    	if isValidRun == true 
			fin_dot_prod_opt 	= String
			fin_spmv_opt 		= String
			fin_mg_opt_nt 		= String 
			fin_waxpby_opt 		= String
			fin_action 			= String 
      		fin_hpcg_validity 	= "HPCG result is VALID with a GFLOP/s rating of $(totalGflops) HPCG 2.4 rating for historical reasons is $(totalGflops24)\n"
			
      		if (A.is_dot_prod_optimized==false) 
        		fin_dot_prod_opt	= "Reference version of ComputeDotProduct used. Performance results are most likely suboptimal\n"
			else
				fin_dot_prod_opt	= "A.is_dot_prod_optimized==true\n"
      		end
      		if (A.is_spmv_optimized ==false) 
        		fin_spmv_opt 	= "Reference version of ComputeSPMV used. Performance results are most likely suboptimal\n"
      		end
      		if (A.is_mg_optimized == false) 
        		if (A.geom.numThreads>1)
          			fin_mg_opt_nt 	= "Reference version of ComputeMG used and number of threads greater than 1 . Performance results are severely suboptimal\n"
        		else # numThreads ==1
          			fin_mg_opt_nt 	= "Reference version of ComputeMG used . Performance results are most likely suboptimal\n"
      			end
			else
				 fin_mg_opt_nt = "A.is_mg_optimized == true\n"
	  		end
      		if (A.is_waxpby_optimized ==false) 
        		fin_waxpby_opt 	= "Reference version of ComputeWAXPBY used. Performance results are most likely suboptimal\n"
			else
				fin_waxpby_opt = "A.is_waxpby_optimized == true\n"
      		end
      		if (times[1]>=minOfficialTime) 
        		fin_action = "Please upload results from the YAML file contents to http:#hpcg-benchmark.org\n"
      		else 
		  		fin_action = "Results are valid but execution time (sec) is $(times[1])\n"
						  #=
        		if (quickPath) 
          			doc.get("Final Summary").add("You have selected the QuickPath option", "Results are official for legacy installed systems with confirmation from the HPCG Benchmark leaders.")
          			doc.get("Final Summary").add("After confirmation please upload results from the YAML file contents to","http:#hpcg-benchmark.org")
        		else 
          			doc.get("Final Summary").add("Official results execution time (sec) must be at least",minOfficialTime)
				end
			=#
	  		end
			println(report_result, final_summary_title, fin_hpcg_validity, fin_dot_prod_opt, fin_spmv_opt, fin_mg_opt_nt, fin_waxpby_opt, fin_action)

     	else 
			
      		fin_hpcg_validity = "HPCG result is","INVALID.\n"
      		fin_action = "Please review the YAML file contents You may NOT submit these results for consideration.\n"
			println(report_result, final_summary_title, fin_hpcg_validity, fin_action)
	  	end
		close(report_result)

   end 
end
