#=
 @file SetupHalo_ref.cpp

 HPCG routine
=#

using MPI
include("SpMatrix.jl")
include("Geometry.jl")
#include("SetupHalo.jl")
#include("mytimer.jl")

#=
  Reference version of SetupHalo that prepares system matrix data structure and creates data necessary
  for communication of boundary values of this process.

  @param[inout] A    The known system matrix

  @see ExchangeHalo
=#
function setup_halo_ref!(A)

  @debug("SetupHalo_ref")

  #Extract Matrix pieces

  localNumberOfRows = A.localNumberOfRows
  nonzerosInRow     = A.nonzerosInRow
  mIG               = A.mtxIndG
  mIL               = A.mtxIndL
  mtxIndG           = permutedims(reshape(hcat(mIG...), (length(mIG[1]), length(mIG))))
  mtxIndL           = permutedims(reshape(hcat(mIL...), (length(mIL[1]), length(mIL))))

  # LNR = ,"localNumberOfRows," dimsMIG = ,",size(mtxIndG)," dimMIL = ",size(mtxIndL),".")
  for i=1:localNumberOfRows 
    cur_nnz = nonzerosInRow[i]
    for j=1:cur_nnz
	mtxIndL[i,j] = mtxIndG[i,j]
    end
  end

  # Scan global IDs of the nonzeros in the matrix.  Determine if the column ID matches a row ID.  If not:
  # 1) We call the ComputeRankOfMatrixRow function, which tells us the rank of the processor owning the row ID.
  #  We need to receive this value of the x vector during the halo exchange.
  # 2) We record our row ID since we know that the other processor will need this value from us, due to symmetry.

  sendList           = Dict()
  receiveList        = Dict()
  externalToLocalMap = Dict()

  #  TODO: With proper critical and atomic regions, this loop could be threaded, but not attempting it at this time
  for i=1:localNumberOfRows 
    currentGlobalRow = AA.localToGlobalMap[i]
    for j=1:nonzerosInRow[i]
      curIndex = mtxIndG[i,j]
      rankIdOfColumnEntry = ComputeRankOfMatrixRow(A.geom, curIndex)
#      @debug("rank, row , col, globalToLocalMap[col] = ",A.geom.rank, currentGlobalRow, curIndex ,AA.globalToLocalMap[curIndex],"\n")
      if A.geom.rank!=rankIdOfColumnEntry #If column index is not a row index, then it comes from another processor
        receiveList[rankIdOfColumnEntry]=curIndex 
        sendList[rankIdOfColumnEntry]=currentGlobalRow # Matrix symmetry means we know the neighbor process wants my value
      end
    end
  end

  #Count number of matrix entries to send and receive
  totalToBeSent = 0
  for (k,v) in sendList  
    totalToBeSent += length(v)
  end
  totalToBeReceived = 0
  for (k,v) in receiveList 
    totalToBeReceived += length(v)
  end

  # These are all attributes that should be true, due to symmetry
#  @debug("totalToBeSent = $totalToBeSent totalToBeReceived = $totalToBeReceived")
  @assert(totalToBeSent==totalToBeReceived) # Number of sent entry should equal number of received
  @assert(length(sendList)==length(receiveList)) # Number of send-to neighbors should equal number of receive-from
  # Each receive-from neighbor should be a send-to neighbor, and send the same number of entries
  for (k,v) in receiveList
    @assert haskey(sendList,k)
    @assert length(sendList[k])==length(receiveList[k])
  end

  #Build the arrays and lists needed by the ExchangeHalo function.
  sendBuffer        = Array{Float64}(undef,totalToBeSent)
  elementsToSend    = Array{Int64}(undef,totalToBeSent)
  neighbors         = Array{Int64}(undef,length(collect(keys(sendList))))
  receiveLength     = Array{Int64}(undef, length(collect(keys(receiveList))))
  sendLength        = Array{Int64}(undef, length(collect(keys(receiveList))))
  neighborCount     = 1
  receiveEntryCount = 1
  sendEntryCount    = 1

  for (k,v) in receiveList 
    neighborId = k #rank of current neighbor we are processing
    neighbors[neighborCount] = neighborId # store rank ID of current neighbor
    receiveLength[neighborCount] = length(receiveList[neighborId])
    sendLength[neighborCount] = length(sendList[neighborId]) # Get count if sends/receives    
    neighborCount +=1
    
    for i in receiveList[neighborId]
      externalToLocalMap[i] = localNumberOfRows + receiveEntryCount # The remote columns are indexed at end of internals
      receiveEntryCount += 1
    end
    for i in sendList[neighborId] 
      elementsToSend[sendEntryCount] = AA.globalToLocalMap[i] # store local ids of entry to send
      sendEntryCount +=1
    end
  end

  # Convert matrix indices to local IDs
  for i=1:localNumberOfRows
    for j=1:nonzerosInRow[i]
      curIndex = mtxIndG[i,j]
      rankIdOfColumnEntry = ComputeRankOfMatrixRow(A.geom, curIndex)
      if A.geom.rank==rankIdOfColumnEntry # My column index, so convert to local index
        mtxIndL[i,j] = AA.globalToLocalMap[curIndex]
      else # If column index is not a row index, then it comes from another processor
        #mtxIndL[i,j] = externalToLocalMap[curIndex]
      end
    end
  end

  # Store contents in our matrix struct
  numberOfExternalValues = length(externalToLocalMap)
  localNumberOfColumns   = AA.localNumberOfRows + numberOfExternalValues
  numberOfSendNeighbors  = length(sendList)
  totalToBeSent          = totalToBeSent
  elementsToSend         = elementsToSend
  neighbors              = neighbors
  receiveLength          = receiveLength
  sendLength             = sendLength
  sendBuffer             = sendBuffer

  A.localNumberOfColumns   = localNumberOfColumns
  A.numberOfExternalValues = numberOfExternalValues
  A.numberOfSendNeighbors  = numberOfSendNeighbors
  A.totalToBeSent          = totalToBeSent
  A.elementsToSend         = elementsToSend
  A.neighbors              = neighbors
  A.receiveLength          = receiveLength
  A.sendLength             = sendLength
  A.sendBuffer             = sendBuffer

  @debug(" For rank $A.geom.rank of $A.geom.size number of neighbors $A.numberOfSendNeighbors")
  for i = 1:numberOfSendNeighbors
    @debug("     rank = ", A.geom.rank," neighbor = ",neighbors[i]," send/recv length = ", sendLength[i]/receiveLength[i],".")
    for j = 1:sendLength[i]
      @debug("       rank = ", A.geom.rank," elementsToSend[$j] =", elementsToSend[j],".")
    end
  end

end
