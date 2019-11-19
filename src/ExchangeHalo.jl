 #=
@file ExchangeHalo.cpp

 HPCG routine
=#
using MPI
# Compile this routine only if running with MPI

include("Geometry.jl")


#=
  Communicates data that is at the border of the part of the domain assigned to this processor.

  @param[in]    A The known system matrix
  @param[inout] x On entry: the local vector entries followed by entries to be communicated on exit: the vector with non-local entries updated by other processors
=#
#fix call
function exchange_halo!(x, A) 

  # Extract Matrix pieces

 if MPI.Initialized()== true
  localNumberOfRows = A.localNumberOfRows
  num_neighbors     = A.numberOfSendNeighbors
  receiveLength     = A.receiveLength
  sendLength        = A.sendLength
  neighbors         = A.neighbors
  sendBuffer        = A.sendBuffer
  totalToBeSent     = A.totalToBeSent
  elementsToSend    = A.elementsToSend

  size = Int64
  rank = Int64 # Number of MPI processes, My process ID

  size = MPI.Comm_size(MPI.COMM_WORLD)
  rank = MPI.Comm_rank(MPI.COMM_WORLD)
  

  #
  #  first post receives, these are immediate receives
  #  Do not wait for result to come, will do that at the
  #  wait call below.
  #

  MPI_MY_TAG = 99
  #array of requests
  request = Array{MPI.Request}(undef, num_neighbors)

  #
  # Externals are at end of locals
  #

  # x_external = x[localNumberOfRows+1:length(x)]

  # In C, x_external points to the (localNumberOfRows)th
  # element of x. The point is to receive indices from neighbors
  # in buff(an array) so that these local indices are stored at the 
  # end of local indices
  
  # In Julia we dont need to preserve this marker for the 
  # 

  # Post receives first
  # TODO: Thread this loop


  if A.geom.rank == 0
	  fs =  open("mpi_ircv_0.txt", "a")
  else
 	  fs =  open("mpi_ircv_1.txt", "a")
  end
  println(fs, "break")



  for i = 1 : num_neighbors
    n_recv 		= receiveLength[i]
    buff 		= Array{Int64,1}(undef,n_recv)

    fill!(buff,0) 			# because buff has double values 

    request[i]  = MPI.Irecv!(buff, neighbors[i], MPI_MY_TAG, MPI.COMM_WORLD)

    x		 	= vcat(x, buff) 			# vector assignment is ok because 
											# x is not a primitive datatype
  end

  println(fs, "$(x[length(x)])  ")

  #
  # Fill up send buffer
  #

  # TODO: Thread this loop
  for i=1:totalToBeSent
	sendBuffer[i] = x[elementsToSend[i]]
  end

  #
  # Send to each neighbor
  #

  # TODO: Thread this loop

  if A.geom.rank == 0
	  fs =  open("mpi_send_0.txt", "a")
  else
 	  fs =  open("mpi_send_1.txt", "a")
  end
  println(fs, "break")

  for i = 1:num_neighbors
    n_send = sendLength[i]
	MPI.Send(sendBuffer[1:n_send], neighbors[i], MPI_MY_TAG, MPI.COMM_WORLD)
	println(fs, "$(sendBuffer[length(sendBuffer)])")
  end



  #
  # Complete the reads issued above
  #

  status = MPI.Status

  #  status = Ref{Status}() as in MPI.jl example package
  # TODO: Thread this loop
  # BUG_INFO 
  # -Check which process was not able to get out of IRecv [ line 61 ]
  # Both process stall at MPI.Wait()
  #
  # Bug occurs in the OptimizedCG Phase 
  # in main.jl at 302 when calling cg!
  #
  # main > cg > mg > symgs > exchangeHalo()
  #
  # -Confirm ExchangeHalo behaviour
  # Exchange Halo is not working as expected
  # all values recieved by every process is 0 if init of buff is fill!(0)(Julia).
  # for process 1 values rcvd other than 0 are negative.
  # most values sent by every process is 0, 1 or NaN (Julia).
  # -exchnageHalo is not working as expected after pointer arithmetic is corrected:
  # HOW
  # 	pointer arithmetic in x_external was creating a copy of the
  # 	[in] vector x in Julia.
  # NOW 
  # x = [local...indices, num_neighbor_1_external_indices, ..., num_beighbor_size-1_external_indices]


  for i = 1: num_neighbors

    if MPI.Wait!(request[i]) == 1
	  println("process number $(A.geom.rank)")
      exit(-1) #TODO: have better error exit
    end

  end
  A.sendBuffer = sendBuffer
  request = nothing
  end
  return
end 
