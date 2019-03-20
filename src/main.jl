#=
 @file main.cpp

 HPCG routine
=#
#   Main routine of a program that calls the HPCG conjugate gradient
#   solver to solve the problem, and then prints results.

import MPI
include("hpcg.jl")

include("CheckAspectRatio.jl")
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
include("SpMatrix.jl")
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
function main() 
  MPI.Init()

  params=HPCG_Params

  HPCG_Init(length(ARGS), ARGS, params)

  #Check if QuickPath option is enabled.
  #If the running time is set to zero, we minimize all paths through the program

  quickPath = (params.runningTime==0)

  size = params.comm_size
  rank = params.comm_rank # Number of MPI processes, My process ID

  if size < 100 && rank==0 
	@debug("Process ",rank, " of ",size," is alive with " , params.numThreads , " threads.") 
  end 

  if rank==0 
    c = readline()
  end

  MPI.Barrier(MPI.COMM_WORLD)

  nx = Int64
  ny = Int64
  nz = Int64
  nx = params.nx
  ny = params.ny
  nz = params.nz
  ierr = 0  #Used to check return codes on function calls

  ierr = CheckAspectRatio(0.125, nx, ny, nz, "local problem", rank==0)

  if ierr!=0
    return ierr
  end

   #########################
   ## Problem setup Phase ##
   #########################

  t1 = mytimer() #INCLUDE CORRECT TIMER 

  #Construct the geometry and linear system
  geom = Geometry
  GenerateGeometry(size, rank, params.numThreads, params.pz, params.zl, params.zu, nx, ny, nz, params.npx, params.npy, params.npz, geom)

  ierr = CheckAspectRatio(0.125, geom.npx, geom.npy, geom.npz, "process grid", rank==0)

  if ierr!=0
    return ierr
  end

  # Use this array for collecting timing information
  times =  zeros(10)
  setup_time = mytimer() #INCLUDE CORRECT TIMER 
  A = SpMatrix
  InitializeSparseMatrix(A, geom)
  # b x and xexact . Array{T,1}
  b =  Vector 
  x =  Vector 
  xexact =  Vector

  GenerateProblem(A, b, x, xexact)
  SetupHalo(A)
  numberOfMgLevels = 4 #Number of levels including first
  curLevelMatrix = A

  for level = 1:numberOfMgLevels
    GenerateCoarseProblem(curLevelMatrix)
    curLevelMatrix = curLevelMatrix.Ac #  Make the just-constructed coarse grid the next level
  end

  setup_time = mytimer() - setup_time #Capture total time of setup
  times[9] = setup_time #Save it for reporting
  curLevelMatrix = A
  curb = b
  curx = x
  curxexact = xexact
  for level = 0 : numberOfMgLevels
     CheckProblem(curLevelMatrix, curb, curx, curxexact)
     curLevelMatrix = curLevelMatrix.Ac # Make the nextcoarse grid the next level
     curb = 0 # No vectors after the top level
     curx = 0
     curxexact = 0
  end


  data = CGData
  InitializeSparseCGData(A, data)



  ####################################
  ## Reference SpMV+MG Timing Phase ##
  ####################################

  # Call Reference SpMV and MG. Compute Optimization time as ratio of times in these routines

  nrow = A.localNumberOfRows
  ncol = A.localNumberOfColumns

  x_overlap = Vector(undef, ncol) #  Overlapped copy of x vector
  b_computed = Vector(undef, nrow) #  Computed RHS vector
 ####################################
#  InitializeVector(x_overlap, ncol)
#  InitializeVector(b_computed, nrow) 
#####################################


  # Record execution time of reference SpMV and MG kernels for reporting times
  # First load vector with random values
  FillRandomVector(x_overlap)

  numberOfCalls = 10
  if quickPath 
	numberOfCalls = 1 #QuickPath means we do on one call of each block of repetitive code
  end
  t_begin = mytimer()
  for i= 1: numberOfCalls
    ierr = ComputeSPMV_ref(A, x_overlap, b_computed) # b_computed = A*x_overlap
    if ierr 
	@debug ("Error in call to SpMV: $ierr .\n")
    end
    ierr = ComputeMG(A, b_computed, x_overlap) # b_computed = Minv*y_overlap
    if (ierr) 
	@debug ("Error in call to MG: $ierr .\n") 
    end
  end 
  times[8] = (mytimer() - t_begin)/( numberOfCalls) # Total time divided by number of calls.
  if rank==0
	 @debug("Total SpMV+MG timing phase execution time in main (sec) = $(mytimer()-t1)\n")
  end 

  ###############################
  ## Reference CG Timing Phase ##
  ###############################

  t1 = mytimer()
  global_failure = 0 # assume all is well: no failures

  niters = 0
  totalNiters_ref = 0
  normr = 0.0
  normr0 = 0.0
  refMaxIters = 50
  numberOfCalls = 1 # Only need to run the residual reduction analysis once

  # Compute the residual reduction for the natural ordering and reference kernels
  ref_times = zeros(9)
  tolerance = 0.0 # Set tolerance to zero to make all runs do maxIters iterations
  err_count = 0
  for i= 1:numberOfCalls
    x = zeros(length(x))
    ierr = CG_ref( A, data, b, x, refMaxIters, tolerance, niters, normr, normr0, ref_times[0], true)
    if ierr
	 err_count+=1 # count the number of errors in CG
    end
    totalNiters_ref += niters
  end
  if rank == 0 && err_count
	 @debug("$err_count error(s) in call(s) to reference CG.")
  end
  refTolerancei:Float64 = normr / normr0

  #Call user-tunable set up function.
  t7 = mytimer()
  OptimizeProblem(A, data, b, x, xexact)
  t7 = mytimer() - t7
  times[7] = t7
  if rank==0
	 @debug("Total problem setup time in main (sec) = $(mytimer()-t1)") 
  end
  if geom.size == 1
	 WriteProblem(geom, A, b, x, xexact)
  end


  ##############################
  ## Validation Testing Phase ##
  ##############################

  t1 = mytimer()
  testcg_data=TestCGData 
  testcg_data.count_pass = testcg_data.count_fail = 0
  TestCG(A, data, b, x, testcg_data)

  testsymmetry_data = TestSymmetryData 
  TestSymmetry(A, b, xexact, testsymmetry_data)

  if (rank==0) 
	@debug("Total validation (TestCG and TestSymmetry) execution time in main (sec) = $(mytimer() - t1)")
  end

  t1 = mytimer()

  ##############################
  ## Optimized CG Setup Phase ##
  ##############################

  niters = 0
  normr = 0.0
  normr0 = 0.0
  err_count = 0
  tolerance_failures = 0

  optMaxIters = 10*refMaxIters
  optNiters = refMaxIters
  opt_worst_time = 0.0

  opt_times =  zeros(9)

  # Compute the residual reduction and residual count for the user ordering and optimized kernels.
  for i=1:numberOfCalls
    x = zeros(length(x)) # start x at all zeros
    last_cummulative_time = opt_times[0]
    ierr = CG( A, data, b, x, optMaxIters, refTolerance, niters, normr, normr0, opt_times[0], true)
    if ierr 
	err_count +=1 # count the number of errors in CG
    end
    if normr / normr0 > refTolerance
	 tolerance_failures+=1 # the number of failures to reduce residual
    end
    # pick the largest number of iterations to guarantee convergence
    if niters > optNiters
	optNiters = niters
    end
    current_time = opt_times[0] - last_cummulative_time
    if current_time > opt_worst_time
	 opt_worst_time = current_time
    end
  end

# Get the absolute worst time across all MPI ranks (time in CG can be different)
  local_opt_worst_time = opt_worst_time
  MPI.Allreduce(local_opt_worst_time, opt_worst_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD)


  if rank == 0 && err_count
	 @debug("$err_count  error(s) in call(s) to optimized CG.") 
  end
  if tolerance_failures 
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

  total_runtime = params.runningTime
  numberOfCgSets = int(total_runtime / opt_worst_time) + 1 # Run at least once, account for rounding

  if rank==0 
    @debug("Projected running time: $total_runtime seconds") 
    @debug("Number of CG sets: $numberOfCgSets") 
  end

  # This is the timed run for a specified amount of time. 

  optMaxIters = optNiters
  optTolerance = 0.0  # Force optMaxIters iterations
  testnorms_data = TestNormsData
  testnorms_data.samples = numberOfCgSets
  testnorms_data.values = Array{Float64}(undef, numberOfCgSets)

  for i=1: numberOfCgSets
    x = zeros(length(x)) # Zero out x
    ierr = CG( A, data, b, x, optMaxIters, optTolerance, niters, normr, normr0, times[0], true)
    if ierr 
	@debug("Error in call to CG: $ierr.\n") 
    end
    if rank==0 
	@debug("Call [$i] Scaled Residual [$(normr/normr0)]") 
    end
    testnorms_data.values[i] = normr/normr0 # Record scaled residual from this run
  end

  # Compute difference between known exact solution and computed solution
  # All processors are needed here.
  residual::Float64 = 0
  ierr = ComputeResidual(A.localNumberOfRows, x, xexact, residual)
  if ierr 
	@debug("Error in call to compute_residual: $ierr.\n") 
  end
  if rank==0 
	@debug("Difference between computed and exact  = $residual .\n") 
  end

  ## Test Norm Results
  ierr = TestNorms(testnorms_data)

  ####################
  ## Report Results ##
  ####################

  # Report results to YAML file
  ReportResults(A, numberOfMgLevels, numberOfCgSets, refMaxIters, optMaxIters, times[0], testcg_data, testsymmetry_data, testnorms_data, global_failure, quickPath)

  #Clean up
  A = nothing # This delete will recursively delete all coarse grid data
  data =  nothing
  x =  nothing
  b = nothing 
  xexact = nothing
  x_overlap = nothing
  b_computed = nothing
  testnorms_data.values = nothing



  HPCG_Finalize()

  # Finish up
  MPI.Finalize()
  return 0
end
main()
