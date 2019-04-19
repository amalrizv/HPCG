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
    mIG               = A.mtxIndG
    mIL               = A.mtxIndL
    mtxIndG           = permutedims(reshape(hcat(mIG...), (length(mIG[1]), length(mIG))))
    mtxIndL           = permutedims(reshape(hcat(mIL...), (length(mIL[1]), length(mIL))))

    # In the non-MPI case we simply copy global indices to local index storage
    @static if MPI.Initialized() == 0
        # LNR = ,"localNumberOfRows," dimsMIG = ,",size(mtxIndG)," dimMIL = ",size(mtxIndL),".")
        for i = 1:localNumberOfRows 
            cur_nnz = nonzerosInRow[i]
            for j = 1:cur_nnz
                mtxIndL[i,j] = mtxIndG[i,j]
            end
        end
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
                #      @debug("rank, row , col, globalToLocalMap[col] = ",A.geom.rank, currentGlobalRow, curIndex ,AA.globalToLocalMap[curIndex],"\n")
                if  A.geom.rank != rankIdOfColumnEntry # If column index is not a row index, then it comes from another processor
                    if haskey(receiveList, rankIdOfColumnEntry)
                        push!(receiveList[rankIdOfColumnEntry], curIndex) 
                    else
                        receiveList[rankIdOfColumnEntry] = Set(curIndex)
                    end

                    if haskey(sendList, rankIdOfColumnEntry)
                        push!(sendList[rankIdOfColumnEntry], curIndex) 
                    else
                        sendList[rankIdOfColumnEntry] = Set(curIndex)
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

        #Build the arrays and lists needed by the ExchangeHalo function.
        sendBuffer        = Array{Float64}(undef, totalToBeSent)
        elementsToSend    = Array{Int64}(undef, totalToBeSent)
        neighbors         = Array{Int64}(undef, length(collect(keys(sendList))))
        receiveLength     = Array{Int64}(undef, length(collect(keys(receiveList))))
        sendLength        = Array{Int64}(undef, length(collect(keys(receiveList))))
        neighborCount     = 1
        receiveEntryCount = 1
        sendEntryCount    = 1

        for (k,v) in receiveList 
            neighborId                   = k # rank of current neighbor we are processing
            neighbors[neighborCount]     = neighborId # store rank ID of current neighbor
            receiveLength[neighborCount] = length(receiveList[neighborId])
            sendLength[neighborCount]    = length(sendList[neighborId]) # Get count if sends/receives

            for i in receiveList[neighborId]
                externalToLocalMap[i] = localNumberOfRows + receiveEntryCount # The remote columns are indexed at end of internals
                receiveEntryCount += 1
            end

            for i in sendList[neighborId] 
                elementsToSend[sendEntryCount] = A.globalToLocalMap[i] # store local ids of entry to send
                sendEntryCount += 1
            end

            neighborCount +=1
        end

        # Convert matrix indices to local IDs
        for i = 1:localNumberOfRows
            for j = 1:nonzerosInRow[i]
                curIndex = mtxIndG[i,j]
                #curIndex =  27 #cannot find this index in A.globalToLocalMap
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
end
