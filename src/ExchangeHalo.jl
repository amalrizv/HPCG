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
  x_external = x[localNumberOfRows:length(x)]

  # Post receives first
  # TODO: Thread this loop
  for i = 1 : num_neighbors
    n_recv = receiveLength[i]
    buff = Array{Int64,1}(undef,n_recv)
    request[i]  = MPI.Irecv!(buff, neighbors[i], MPI_MY_TAG, MPI.COMM_WORLD)
    vcat(x_external, buff)
  end


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
  #AMAL:TODO: fix pointer arithmetic 
  for i = 1:num_neighbors
    n_send = sendLength[i]
    start = (i-1)n_send+1
    finish = start+n_send-1
    MPI.Send(sendBuffer[start:finish], neighbors[i], MPI_MY_TAG, MPI.COMM_WORLD)
    
  end

  #
  # Complete the reads issued above
  #

  status = MPI.Status
  # TODO: Thread this loop
  for i = 1: num_neighbors
    if MPI.Wait!(request[i]) == 1
      exit(-1) #TODO: have better error exit
    end
  end

  request = nothing
  end
  return
end
