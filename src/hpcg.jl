#=
  @file hpcg.hpp

  HPCG data structures and functions
=#

include("Geometry.jl")


mutable struct HPCG_Params
	 comm_size::Int64 	#Number of MPI processes in MPI_COMM_WORLD
	 comm_rank::Int64 	#This process' MPI rank in the range [0 to comm_size - 1]
	 numThreads::Int64  	#This process' number of threads
	 nx::Int64  		#Number of processes in x-direction of 3D process grid
	 ny::Int64 		#Number of processes in y-direction of 3D process grid
	 nz::Int64 		#Number of processes in z-direction of 3D process grid
	 runningTime::Int64  	#Number of seconds to run the timed portion of the benchmark
	 npx::Int64  		#Number ::Int64 of x-direction grid pos for each local subdomain
	 npy::Int64 		# Number of y-direction grid pos for each local subdomain
	 npz::Int64  		# Number of z-direction grid pos for each local subdomain
	 pz::Int64  		# Partition in the z processor dimension, default is npz
	 zl::Int64  		#nz for processors in the z dimension with value less than pz
	 zu::Int64  		#nz for processors in the z dimension with value greater than pz
end

