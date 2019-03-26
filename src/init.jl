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

function HPCG_Init(arg_hpcg, params) 
  #arg_hpcg  = [np, nx,ny,nz, d,rt, mpi]
  opts = collect(keys(arg_hpcg))
  vals = collect(values(arg_hpcg))
  #Documentation asks us to only provide nx ny nz np and rt not pz zl zu npx npy npz
  # If these are not provided in command line then assign them zero, therefore last 6 values are 0
  #TODO :Add option for pz zl zu npz npy npz
  #     :Add option to read from File

  iparams = [arg_hpcg["nx"], arg_hpcg["ny"], arg_hpcg["nz"], arg_hpcg["rt"], 0,0,0,0,0,0]
  cparams = ["--nx=", "--ny=", "--nz=", "--rt=", "--pz=", "--zl=", "--zu=", "--npx=", "--npy=", "--npz="]
  nparams = length(cparams)
  broadcastParams = false # Make true if parameters read from file.


  

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
  if arg_hpcg["USE_MPI"] == true
    if broadcastParams == true 
      MPI.Bcast( iparams, nparams, 0, MPI.COMM_WORLD )
    end
  end
  if arg_hpcg["USE_MPI"] == true
    params = HPCG_Params(MPI.size(MPI.COMM_WORLD), MPI.rank(MPI.COMM_WORLD), 1, iparams[1], iparams[2], iparams[3], iparams[4], iparams[8], iparams[9], iparams[10], iparams[5], iparams[6], iparams[7])
  else
    params = HPCG_Params(1, 0, 1, iparams[1], iparams[2], iparams[3], iparams[4], iparams[8], iparams[9], iparams[10], iparams[5], iparams[6], iparams[7])
  end 

  #=
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
  =#


  date = now()
  fname =  "hpcg"*"_"*string(date)*".txt"
  io = open(fname, "w+")
  HPCG_Fout = SimpleLogger(io,Logging.Debug)
  HPCG_fout = global_logger(HPCG_Fout)

  iparams  = nothing

  return params
end
