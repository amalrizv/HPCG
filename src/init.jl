using MPI
using Dates
using Logging

include("hpcg.jl")

include("ReadHpcgDat.jl")

#=Julia has its own function for this
static int
startswith(const char * s, const char * prefix) {
  size_t n = strlen( prefix )
  if (strncmp( s, prefix, n ))
    return 0
  return 1
}
=#
#=
  Initializes an HPCG run by obtaining problem parameters (from a file or
  command line) and then broadcasts them to all nodes. It also initializes
  login I/O streams that are used throughout the HPCG run. Only MPI rank 0
  performs I/O operations.

  The function assumes that MPI has already been initialized for MPI runs.

  @param[in] argc_p the pointer to the "argc" parameter passed to the main() function
  @param[in] argv_p the pointer to the "argv" parameter passed to the main() function
  @param[out] params the reference to the data structures that is filled the basic parameters of the run

  @return returns 0 upon success and non-zero otherwise

  @see HPCG_Finalize
=#

function HPCG_Init(length, ARGS , params) 
  argc = length(ARGS)
  argv = ARGS
  fname = String
  i =Int64 
  j =Int64
  cparams[:,7] = ["--nx=", "--ny=", "--nz=", "--rt=", "--pz=", "--zl=", "--zu=", "--npx=", "--npy=", "--npz="]
  nparams = sizeof(cparams) / sizeof(cparams[0])
  broadcastParams = false # Make true if parameters read from file.

  iparams = Array{Int64,1}(undef, nparams)

  
  # for sequential and some MPI implementations it's OK to read first three args */
  for i = 1:4
    if  ARGS[i]< 10
	 iparams = Int64(ARGS[1])
	 iparams[i] = 0
    end
  end
  # for some MPI environments, command line arguments may get complicated so we need a prefix */
  i = 1 
  while i <= length[ARGS] 
    for j =1:nparams 
      if startswith(ARGS[i], cparams[j])
	#startswith(s::AbstractString, prefix::AbstractString)
	#Returns true if s starts with prefix. If prefix is a vector or set of characters, tests whether the first character of s belongs to that set.
          iparams[j] = 0
      end
    end 
   i = i+1
  end
  # Check if --rt was specified on the command line
  rt  = iparams+3  #Assume runtime was not specified and will be read from the hpcg.dat file
  if ! iparams[3]
	 rt = 0 #If --rt was specified, we already have the runtime, so don't read it from file
  end
  if ! iparams[0] && ! iparams[1] && ! iparams[2] # no geometry arguments on the command line 
    ReadHpcgDat(iparams, rt, iparams+7)
    broadcastParams = true
  end

  #Check for small or unspecified nx, ny, nz values
  #If any dimension is less than 16, make it the max over the other two dimensions, or 16, whichever is largest
  for i =1:3 
    if iparams[i] < 16
      for j = 1:2 
        if iparams[(i+j)%3] > iparams[i]
          iparams[i] = iparams[(i+j)%3]
	end
      end
    end 		
    if iparams[i] < 16
      iparams[i] = 16
    end
  end

#Broadcast values of iparams to all MPI processes
#ifndef HPCG_NO_MPI
  if (broadcastParams) 
    MPI.Bcast( iparams, nparams, MPI_INT, 0, MPI.COMM_WORLD )
  end

  params.nx = iparams[1]
  params.ny = iparams[2]
  params.nz = iparams[3]

  params.runningTime = iparams[4]
  params.pz = iparams[5]
  params.zl = iparams[6]
  params.zu = iparams[7]

  params.npx = iparams[8]
  params.npy = iparams[9]
  params.npz = iparams[10]

  params.comm_rank = 0
  params.comm_size = 1

  params.numThreads = 1
#  params.numThreads = omp_get_num_threads()
  date = Dates.DateTimenow()
  fname =  "hpcg"*"_"*string(date)*".txt"
  io = open(fname, "w+")
  HPCG_fout = SimpleLogger(io)
  Logging.configure(HPCG_fout, level=DEBUG)

  if 0 == params.comm_rank
    Logging.configure(HPCG_fout, filename=fname)
  else 
  fname =  "hpcg"*"_"*string(date)*".txt"
    Logging.configure(HPCG_fout, filename=fname)
  end

  iparams  = nothing

  return 0
end
