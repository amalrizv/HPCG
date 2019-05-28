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
function exchange_halo(A, x) 

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

  xv   = x

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
  reqs = Array{MPI.Request}(undef, num_neighbors)

  #
  # Externals are at end of locals
  #
  x_external = resize!(xv , length(xv)+ localNumberOfRows)

  # Post receives first
  # TODO: Thread this loop
  for i = 1 : num_neighbors
    n_recv = receiveLength[i]
    reqs[i]  = MPI.Irecv!(x_external, neighbors[i], MPI.MY_TAG, MPI.COMM_WORLD)
    x_external += n_recv
  end


  #
  # Fill up send buffer
  #

  # TODO: Thread this loop
  for i=1:totalToBeSent
	sendBuffer[i] = xv[elementsToSend[i]]
  end

  #
  # Send to each neighbor
  #

  # TODO: Thread this loop
  for i = 1:num_neighbors
    n_send = sendLength[i]
    MPI.Send(sendBuffer, neighbours[i], MPI.MY_TAG, MPI.COMM_WORLD)
    sendBuffer += n_send
  end

  #
  # Complete the reads issued above
  #

  status = MPI.Status
  # TODO: Thread this loop
  for i = 1: num_neighbors
    if MPI.Wait(request[i]) 
      exit(-1) #TODO: have better error exit
    end
  end

  request = nothing
  end
  return
end
