 #=
@file ExchangeHalo.cpp

 HPCG routine
=#
# Compile this routine only if running with MPI





#=
  Communicates data that is at the border of the part of the domain assigned to this processor.

  @param[in]    A The known system matrix
  @param[inout] x On entry: the local vector entries followed by entries to be communicated on exit: the vector with non-local entries updated by other processors
=#
#fix call
function exchange_halo!(x::Array{Float64,1} , A::HPCGSparseMatrix) 
  # Extract Matrix pieces
 if MPI.Initialized()== true
  localNumberOfRows = A.localNumberOfRows
  num_neighbors     = A.numberOfSendNeighbors
  receiveLength     = A.receiveLength
  sendLength        = A.sendLength
  neighbors         = A.neighbors
  sendBuffer        = A.sendBuffer
  totalToBeSent     = A.totalToBeSent
  elementsToSend    = A.elementsToSend # corresponds with C version

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

  # Post receives first
  # TODO: Thread this loop

  # any changes in xv should not be reflected in x
  
  offset_start = localNumberOfRows+1  
  for i = 1 : num_neighbors
    n_recv 		= receiveLength[i]
    buff 		= Array{Float64,1}(undef,n_recv)
    offset_stop = offset_start +n_recv-1
    zero_fill!(buff) 			# because buff has double values 

    request[i]  = MPI.Irecv!(buff, neighbors[i], MPI_MY_TAG, MPI.COMM_WORLD)
			# MPI.Irecv! shows  buff array does not receive anything
	
	x[offset_start:offset_stop] = buff
	offset_start = offset_stop+1


  end
# x has both local and external indices.

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

  # Only local indices from x are used here
  for i = 1:num_neighbors
    n_send = sendLength[i]
	MPI.Send(sendBuffer[1:n_send], neighbors[i], MPI_MY_TAG, MPI.COMM_WORLD)
  end




  #
  # Complete the reads issued above
  #

  status = MPI.Status

  #  status = Ref{Status}() as in MPI.jl example package
  # TODO: Thread this loop
  for i = 1: num_neighbors

    if MPI.Wait!(request[i]) == 1
      exit(-1) #TODO: have better error exit
    end

  end
  A.sendBuffer = sendBuffer
  request = nothing
  end
  return
end 
