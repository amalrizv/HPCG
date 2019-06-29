#=
 @file SetupHalo_ref.jl

 HPCG routine
=#

using MPI
include("Geometry.jl")

#=
  Reference version of SetupHalo that prepares system matrix data structure and creates data necessary
  for communication of boundary values of this process.

  @param[inout] A    The known system matrix

  @see ExchangeHalo
=#
function setup_halo_ref!(A)

    @debug("In SetupHalo_ref")

    # Extract Matrix pieces

    localNumberOfRows = A.localNumberOfRows
    nonzerosInRow     = A.nonzerosInRow
    mtxIndG               = A.mtxIndG
    mtxIndL               = A.mtxIndL
    # In the non-MPI case we simply copy global indices to local index storage
    if MPI.Initialized() == false
        # LNR = ,"localNumberOfRows," dimsMIG = ,",size(mtxIndG)," dimMIL = ",size(mtxIndL),".")
        A.mtxIndL = mtxIndG
    else

        # Scan global IDs of the nonzeros in the matrix.  Determine if the column ID matches a row ID.  If not:
        # 1) We call the compute_rank_of_matrix_row function, which tells us the rank of the processor owning the row ID.
        #  We need to receive this value of the x vector during the halo exchange.
        # 2) We record our row ID since we know that the other processor will need this value from us, due to symmetry.

        sendList           = Dict{Int64, Set{Int64}}()
        receiveList        = Dict{Int64, Set{Int64}}()
        externalToLocalMap = Dict{Int64, Int64}()

        #  TODO: With proper critical and atomic regions, this loop could be threaded, but not attempting it at this time
        for i = 1:localNumberOfRows 
            currentGlobalRow = A.localToGlobalMap[i]
            for j = 1:nonzerosInRow[i]
                curIndex            = mtxIndG[i,j]
                rankIdOfColumnEntry = compute_rank_of_matrix_row(A.geom, curIndex)
		@show rankIdOfColumnEntry
                #      @debug("rank, row , col, globalToLocalMap[col] = ",A.geom.rank, currentGlobalRow, curIndex ,AA.globalToLocalMap[curIndex],"\n")
                if  A.geom.rank != rankIdOfColumnEntry # If column index is not a row index, then it comes from another processor

		    if haskey(receiveList, rankIdOfColumnEntry)
		    	push!(receiveList[rankIdOfColumnEntry], curIndex)
 		    else
			receiveList[rankIdOfColumnEntry] = Set(curIndex)
		    end
                    if haskey(sendList, rankIdOfColumnEntry) 
                    	push!(sendList[rankIdOfColumnEntry], currentGlobalRow) 
		    else
                    	sendList[rankIdOfColumnEntry] = Set(currentGlobalRow)
		    end

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
	@show totalToBeSent
	@show totalToBeReceived
        # TODO KCH: the following should only execute if debugging is enabled!
        # These are all attributes that should be true, due to symmetry
        #  @debug("totalToBeSent = $totalToBeSent totalToBeReceived = $totalToBeReceived")
        @assert(totalToBeSent==totalToBeReceived) # Number of sent entry should equal number of received
        @assert(length(sendList)==length(receiveList)) # Number of send-to neighbors should equal number of receive-from
        # Each receive-from neighbor should be a send-to neighbor, and send the same number of entries
        for (k,v) in receiveList
            @assert haskey(sendList,k)
            @assert length(sendList[k])==length(receiveList[k])
        end
        len_snd_list  = length(collect(keys(sendList)))
        len_rcv_list  = length(collect(keys(receiveList)))
        #Build the arrays and lists needed by the ExchangeHalo function.
        sendBuffer        = Array{Float64}(undef, totalToBeSent)
        elementsToSend    = Array{Int64}(undef, totalToBeSent)
        neighbors         = Array{Int64}(undef, len_snd_list)
        receiveLength     = Array{Int64}(undef, len_rcv_list)
        sendLength        = Array{Int64}(undef, len_rcv_list)
        neighborCount     = 1
        receiveEntryCount = 1
        sendEntryCount    = 1

        for (k,v) in receiveList 
            neighborId                   = k # rank of current neighbor we are processing
            neighbors[neighborCount]     = neighborId # store rank ID of current neighbor
            receiveLength[neighborCount] = length(receiveList[neighborId])
            sendLength[neighborCount]    = length(sendList[neighborId]) # Get count if sends/receives
	    n_rcv_id = receiveList[neighborId] 
	    n_snd_id = sendList[neighborId]
	   @show n_rcv_id
            for i in n_rcv_id
                externalToLocalMap[i] = localNumberOfRows + receiveEntryCount # The remote columns are indexed at end of internals
                receiveEntryCount += 1
            end

            for i in n_snd_id
                elementsToSend[sendEntryCount] = A.globalToLocalMap[i] # store local ids of entry to send
                sendEntryCount += 1
            end

            neighborCount +=1
        end

        # Convert matrix indices to local IDs
        for i = 1:localNumberOfRows
            for j = 1:nonzerosInRow[i]
                curIndex = mtxIndG[i,j]
                rankIdOfColumnEntry = compute_rank_of_matrix_row(A.geom, curIndex)
                if A.geom.rank == rankIdOfColumnEntry # My column index, so convert to local index
                    mtxIndL[i,j] = A.globalToLocalMap[curIndex]
                else # If column index is not a row index, then it comes from another processor
                    mtxIndL[i,j] = externalToLocalMap[curIndex]
                end
            end
       end

        # Store contents in our matrix struct
        numberOfExternalValues = length(externalToLocalMap)
        localNumberOfColumns   = A.localNumberOfRows + numberOfExternalValues
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

    end # ! NO_MPI
    println("Out of SetupHalo_ref")
    return A
end
