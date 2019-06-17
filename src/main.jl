#=
 @file main.jl

 HPCG routine
=#

#   Main routine of a program that calls the HPCG conjugate gradient
#   solver to solve the problem, and then prints results.

# KCH TODO: only do this if mpi is enabled!
using MPI
include("hpcg.jl")
include("init.jl")
include("CheckAspectRatio.jl")
include("SpMatrix.jl")
include("GenerateGeometry.jl")
include("GenerateProblem.jl")
include("GenerateCoarseProblem.jl")
include("SetupHalo.jl")
include("CheckProblem.jl")
include("ExchangeHalo.jl")
include("OptimizeProblem.jl")
#include("WriteProblem.jl")
#include("ReportResults.jl")
#include("mytimer.jl")
include("ComputeSPMV_ref.jl")
include("ComputeMG.jl")
include("ComputeResidual.jl")
include("CG.jl")
include("CG_ref.jl")
include("Geometry.jl")
#include("Vector.jl")
include("CGData.jl")
include("TestCG.jl")
include("TestSymmetry.jl")
include("TestNorms.jl")

#=
  Main driver program: Construct synthetic problem, run V&V tests, compute benchmark parameters, run benchmark, report results.

  @param[in]  argc Standard argument count.  Should equal 1 (no arguments passed in) or 4 (nx, ny, nz passed in)
  @param[in]  argv Standard argument array.  If argc==1, argv is unused.  If argc==4, argv[1], argv[2], argv[3] will be interpreted as nx, ny, nz, resp.

  @return Returns zero on success and a non-zero value otherwise.

=#



function main(hpcg_args) 

    # retrieve arguments:
    opts = collect(keys(hpcg_args))
    vals = collect(values(hpcg_args))

    if hpcg_args["np"] > 1
        MPI.Init()
	println("initialized")
        hpcg_args["np"] = MPI.Comm_size(MPI.COMM_WORLD)
    end

    params = HPCG_Params
    params = hpcg_init!(params, hpcg_args)

    # Check if QuickPath option is enabled.
    # If the running time is set to zero, we minimize all paths through the program
    quickPath::Bool = (params.runningTime == 0)

    size = params.comm_size
    rank = params.comm_rank # Number of MPI processes, My process ID

    if size < 100 && rank == 0
        @debug("Process $rank of size $size is alive with $(params.numThreads) threads") 
    end

    if MPI.Initialized()
        MPI.Barrier(MPI.COMM_WORLD)
    end

    nx = params.nx
    ny = params.ny
    nz = params.nz

    ierr = 0  # Used to check return codes on function calls
    ierr = check_aspect_ratio(0.125, nx, ny, nz, "local problem", rank == 0)

    if ierr != 0
        return ierr
    end

    #########################
    ## Problem Setup Phase ##
    #########################

    t1 = time_ns() # TODO: INCLUDE CORRECT TIMER 

    # Construct the geometry and linear system
    geom = generate_geometry!(size, rank, params.numThreads, params.pz, params.zl, params.zu, nx, ny, nz, params.npx, params.npy, params.npz)

    ierr = check_aspect_ratio(0.125, geom.npx, geom.npy, geom.npz, "process grid", rank==0)

    if ierr != 0
        return ierr
    end

    # Use this array for collecting timing information
    times      = zeros(10)
    setup_time = time_ns() # TODO: INCLUDE CORRECT TIMER
    A          = initialize_sparse_matrix(geom)

    b, x, xexact = generate_problem!(A)	
    A   = setup_halo!(A)	
    num_mg_levels  = 4 #Number of levels including first

    cur_level_matrix::HPCGSparseMatrix = A
    for level = 1:num_mg_levels-1
        @debug("Generating course problem for level=$level")
        cur_level_matrix  = generate_coarse_problem!(cur_level_matrix) 	
        cur_level_matrix = cur_level_matrix.Ac 		#  Make the just-constructed coarse grid the next level
    end
    @debug("All levels generated")

    setup_time = time_ns() - setup_time #Capture total time of setup
    # TODO: Why is the below commented out?
    #  times[9] = setup_time #Save it for reporting

    cur_level_matrix = A
    curb             = b
    curx             = x
    curxexact        = xexact

    for level = 1:num_mg_levels-1

        check_problem(cur_level_matrix, curb, curx, curxexact)
        cur_level_matrix = A.Ac # Make the nextcoarse grid the next level
        curb           = 0    # No vectors after the top level
        curx           = 0
        curxexact      = 0
    end

    data = CGData
    data = InitializeSparseCGData(A, data)

    ####################################
    ## Reference SpMV+MG Timing Phase ##
    ####################################

    # Call Reference SpMV and MG. Compute Optimization time as ratio of times in these routines

    nrow = A.localNumberOfRows
    ncol = A.localNumberOfColumns
    x_overlap  = Vector{Float64}(undef, ncol) #  Overlapped copy of x vector
    b_computed = Vector{Float64}(undef, nrow) #  Computed RHS vector

    # Record execution time of reference SpMV and MG kernels for reporting times
    # First load vector with random values
#    x_overlap = fill_random_vector!(x_overlap)
     fill!(x_overlap,1)
    num_calls = 10
    if quickPath ==1 
        num_calls = 1 # QuickPath means we do on one call of each block of repetitive code
    end

    t_begin = time_ns()

    for i = 1:num_calls
        ierr, b_computed = compute_spmv_ref!(b_computed,A, x_overlap) # b_computed = A*x_overlap
        #wrong 
	@show b_computed[nrow]
        if ierr != 0
            @error("Error in call to SpMV: $ierr .\n")
        end

        ierr, x_overlap= compute_mg!(x_overlap, A, b_computed) # b_computed = Minv*y_overlap

        if ierr != 0
            @error("Error in call to MG: $ierr .\n") 
        end
    end 

    times[8] = (time_ns() - t_begin)/num_calls # Total time divided by number of calls.

    #@debug begin
    #    if rank == 0
    #        @debug("Total SpMV+MG timing phase execution time in main (sec) = $(time_ns()-t1)\n")
    #    end 
    #end

    ###############################
    ## Reference CG Timing Phase ##
    ###############################

    t1 = time_ns()

    global_failure  = 0 # assume all is well: no failures
    niters          = 0
    totalNiters_ref = 0
    normr           = 0.0
    normr0          = 0.0
    refMaxIters     = 50
    num_calls       = 1 # Only need to run the residual reduction analysis once

    # Compute the residual reduction for the natural ordering and reference kernels
    ref_times = zeros(9)
    tolerance = 0.0 # Set tolerance to zero to make all runs do maxIters iterations
    err_count = 0

    for i = 1:num_calls
        x = zeros(length(x))
        @debug("In     ## Reference CG Timing Phase ## ")
        ierr, A, data, x, niters, normr, normr0, ref_add = cg_ref!(A, data, b, x, refMaxIters, tolerance, niters, normr, normr0, ref_times, true)
        ref_times = ref_add + ref_times

        if ierr != 0
            err_count += 1 # count the number of errors in CG
        end

        totalNiters_ref += niters
    end

    #@debug begin
    #    if rank == 0 
    #        @debug("$err_count error(s) in call(s) to reference CG.")
    #    end
    #end

    refTolerance::Float64 = normr / normr0

    # Call user-tunable set up function.
    t7 = time_ns()
    optimize_problem(A, data, b, x, xexact)
    t7 = time_ns() - t7
    times[7] = t7

    #@debug begin
    #    if rank==0
    #        @debug("Total problem setup time in main (sec) = $(time_ns()-t1)") 
    #    end
    #end

    if geom.size == 1
        # WriteProblem(geom, A, b, x, xexact)
    end

    ##############################
    ## Validation Testing Phase ##
    ##############################

    t1 = time_ns()

    count_pass = 0
    count_fail = 0

    @debug "Length of solution vector: " length(x)

    TestCGdata, times = test_cg!(A, data, b, x, count_pass, count_fail)

    testsymmetry_data = TestSymmetryData 
#    test_symmetry(A, b, xexact, testsymmetry_data)

    if (rank==0) 
        @debug "Total validation (TestCG and TestSymmetry) execution time in main (sec) = " (time_ns() - t1)
    end

    t1 = time_ns()

    ##############################
    ## Optimized CG Setup Phase ##
    ##############################

    niters             = 0
    normr              = 0.0
    normr0             = 0.0
    err_count          = 0
    tolerance_failures = 0

    optMaxIters    = 10*refMaxIters
    optNiters      = refMaxIters
    opt_worst_time = 0.0

   opt_times = zeros(9)

    # Compute the residual reduction and residual count for the user ordering and optimized kernels.
    for i=1:num_calls

        x = zeros(length(x)) # start x at all zeros
        last_cummulative_time = opt_times[1]

        ierr, A, data, x, niters, normr, normr0, opt_add = cg!(A, data, b, x, optMaxIters, refTolerance, niters, normr, normr0, opt_times, true)

        opt_times += opt_add

        if ierr != 0
            err_count +=1 # count the number of errors in CG
        end

        if normr / normr0 > refTolerance
            tolerance_failures += 1 # the number of failures to reduce residual
        end

        # pick the largest number of iterations to guarantee convergence
        if niters > optNiters
            optNiters = niters
        end

        current_time = opt_times[1] - last_cummulative_time

        if current_time > opt_worst_time
            opt_worst_time = current_time
        end
    end

    # Get the absolute worst time across all MPI ranks (time in CG can be different)
    local_opt_worst_time = opt_worst_time
     if MPI.Initialized()==true
	    opt_worst_time = MPI.Allreduce(local_opt_worst_time, MPI.MAX, MPI.COMM_WORLD)
    end

    if rank == 0 && err_count == 1
        @error("$err_count  error(s) in call(s) to optimized CG.") 
    end

    if tolerance_failures == 1 
        global_failure = 1
        if rank == 0
            @debug("Failed to reduce the residual $tolerance_failures times.")
        end
    end

    ###############################
    ## Optimized CG Timing Phase ##
    ###############################

    ## Here we finally run the benchmark phase
    ## The variable total_runtime is the target benchmark execution time in seconds

    total_runtime  = params.runningTime
    numberOfCgSets = floor(total_runtime / opt_worst_time) + 1 # Run at least once, account for rounding
    
    if rank == 0 
        @debug("Projected running time: $total_runtime seconds") 
        @debug("Number of CG sets: $numberOfCgSets") 
    end

    # This is the timed run for a specified amount of time. 

    optMaxIters    = optNiters
    optTolerance   = 0.0  # Force optMaxIters iterations
    vals           = Array{Float64, 1}(undef, Int(numberOfCgSets))
    testnorms_data = TestNormsData(vals, 0.0, 0.0, numberOfCgSets, true)

    for i=1: Int(numberOfCgSets)

        x = zeros(length(x)) # Zero out x
        ierr, A, data, x, niters, normr, normr0, times_add = cg!(A, data, b, x, optMaxIters, optTolerance, niters, normr, normr0, times, true)
        times = times_add+times

        if ierr != 0
            @error("Error in call to CG: $ierr.\n") 
        end

        if rank == 0 
            @debug("Call [$i] Scaled Residual [$(normr/normr0)]") 
        end

        testnorms_data.values[i] = normr/normr0 # Record scaled residual from this run
    end

    # Compute difference between known exact solution and computed solution
    # All processors are needed here.
    residual::Float64 = 0.0
    ierr , residual   = compute_residual!(residual ,A.localNumberOfRows, x, xexact)

    if ierr != 0
        @error("Error in call to compute_residual: $ierr.\n") 
    end

    if rank == 0 
        @debug("Difference between computed and exact  = $residual.\n") 
    end

    ## Test Norm Results
    ierr = test_norms(testnorms_data)

    ####################
    ## Report Results ##
    ####################
    @show length(x)
    # Report results to YAML file
#    times  = ReportResults(A, numberOfMgLevels, numberOfCgSets, refMaxIters, optMaxIters, times, testcg_data, testsymmetry_data, testnorms_data, global_failure, quickPath)
    # Clean up
    A              = nothing # This delete will recursively delete all coarse grid data
    data           = nothing
    x              = nothing
    b              = nothing
    xexact         = nothing
    x_overlap      = nothing
    b_computed     = nothing
    testnorms_data = nothing

    # Finish up
    if MPI.Initialized()==true
	println("message")
  	MPI.Finalize()
    end

    return 0
end
