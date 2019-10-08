#=
 @file SetupHalo_ref.jl

 HPCG routine
=#

using MPI
using DelimitedFiles
include("Geometry.jl")
include("SetupHalo_struct.jl")
#=
  Reference version of SetupHalo that prepares system matrix data structure and creates data necessary
  for communication of boundary values of this process.

  @param[inout] A    The known system matrix

  @see ExchangeHalo
=#
function setup_halo_ref!(A)

#	DEBUG : Same values forwarded from GenerateProblem_ref	
	if A.geom.rank == 1 
		open("mtx_setup_1.txt", "a") do f 
			println(f, A.mtxIndG, A.mtxIndL)
		end
	end
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

		sendList           = Dict{Int64, alt_set}()
		receiveList        = Dict{Int64, alt_set}()
        externalToLocalMap = Dict{Int64, Int64}()
        #  TODO: With proper critical and atomic regions, this loop could be threaded, but not attempting it at this time
        for i = 1:localNumberOfRows 
            currentGlobalRow = A.localToGlobalMap[i] 
            for j = 1:nonzerosInRow[i]
				curIndex            = mtxIndG[i,j]
                rankIdOfColumnEntry = compute_rank_of_matrix_row(A.geom, curIndex)
				for z =0:A.geom.size-1 #create alt_sets for send and rcv lists for all ranks 
					receiveList[z]	= alt_set(0)
					sendList[z]		= alt_set(0)
				end
				## DEBUG
#				open("col_entry_rank_output.txt", "a") do f
#			     	println(f, "rank, row , col, globalToLocalMap[col] = $(A.geom.rank), $currentGlobalRow, $curIndex, $(A.globalToLocalMap[curIndex])")	
#				end
				## DEBUG
				
#					@show A.geom.rank, currentGlobalRow, curIndex, length(A.globalToLocalMap)



            	if  A.geom.rank != rankIdOfColumnEntry # If column index is not a row index, then it comes from another processor
					@show curIndex
					#add to alt_set_add!
					alt_set_add!(receiveList[rankIdOfColumnEntry], curIndex)
					alt_set_add!(sendList[rankIdOfColumnEntry], currentGlobalRow)
					#=
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
					=#
            	end
        	end
    	end
		@show receiveList[0], receiveList[1], sendList[0], sendList[1]
    #Count number of matrix entries to send and receive
    	totalToBeSent = 0
		for (s,t) in sendList
			# sendList is a Dict of Int(s) and alt_set(t)
        	for (k,v) in t.alt
				# alt_set(t)  is a struct of Int(k) and Dict(v)
				totalToBeSent += length(v)
			end
        end

	    totalToBeReceived = 0
		for (s,t) in receiveList
			# receiveList is a Dict of Int(s) and alt_set(t)
	        for (k,v) in t.alt
				# alt_set(t)  is a Dict of Int(k) and Dict(v)
				totalToBeReceived += length(v)
			end
        end
        # TODO KCH: the following should only execute if debugging is enabled!
        # These are all attributes that should be true, due to symmetry
        # @debug("totalToBeSent = $totalToBeSent totalToBeReceived = $totalToBeReceived")
        @assert(totalToBeSent==totalToBeReceived) # Number of sent entry should equal number of received
        @assert(length(sendList)==length(receiveList)) # Number of send-to neighbors should equal number of receive-from
        # Each receive-from neighbor should be a send-to neighbor, and send the same number of entries
        for (k,v) in receiveList
            @assert haskey(sendList,k)
            @assert length(sendList[k].alt)==length(receiveList[k].alt)
        end
        len_snd_list  = length(collect(keys(sendList)))
        len_rcv_list  = length(collect(keys(receiveList)))
        #Build the arrays and lists needed by the ExchangeHalo function.
        sendBuffer        = Array{Float64}(undef, totalToBeSent)
        elementsToSend    = Array{Int64}(undef, totalToBeSent)
        neighbors         = Array{Int64}(undef, len_snd_list)
        receiveLength     = Array{Int64}(undef, len_rcv_list)
        sendLength        = Array{Int64}(undef, len_rcv_list)
        neighborCount     = 0

        for (k,v) in receiveList 
            neighborId                   = k # rank of current neighbor we are processing
            neighbors[neighborCount+1]     = neighborId # store rank ID of current neighbor
            receiveLength[neighborCount+1] = length(receiveList[neighborId].alt)
            sendLength[neighborCount+1]    = length(sendList[neighborId].alt) # Get count if sends/receives
			n_rcv_id = receiveList[neighborId].alt 	# n_rcv_id is a dict of an int and 
													# an alt_set(with fields size(int) 
													# and alt(dict))

	    	n_snd_id = sendList[neighborId].alt		# n_rcv_id is a dict of an int and 
													# an alt_set(with fields size(int)
													# and alt(dict))
@show n_rcv_id
			for(x,y) in n_rcv_id # is a dict of index(key) and factor(value)     		 
				# x(key) is curIndex value 
				# y(value) is the factor to be added
				# Check for i == 16, 48,528 |receiveEntryCount ,0,1,2
		#		if i == 48
		#		  @show localNumberOfRows, 
		#		   exit()
		#	    end
				@show x, localNumberOfRows, y				
				externalToLocalMap[x] = localNumberOfRows + y # The remote columns are indexed at end of internals
            end

			for (x,y) in n_snd_id
				# x(key) is currentGlobalRow value 
				# y(value) is the indexing 

#				if A.geom.rank == 0 
#					open("els_2_send_rank_1.txt", "a") do f 
#						println(f,"$(A.globalToLocalMap[i]) ,$sendEntryCount") 
#					end
#				else
#					open("els_2_send_rank_2.txt", "a") do f 
#						println(f,"$(A.globalToLocalMap[i]) ,$sendEntryCount") 
#					end
#				end

                elementsToSend[y] = A.globalToLocalMap[x] # store local ids of entry to send
            end

            neighborCount +=1
        end

        # Convert matrix indices to local IDs
        for i = 1:localNumberOfRows
            for j = 1:nonzerosInRow[i]
                curIndex = mtxIndG[i,j]
                rankIdOfColumnEntry = compute_rank_of_matrix_row(A.geom, curIndex)
                if A.geom.rank == rankIdOfColumnEntry # My column index, so convert to local index

					# PROBLEM :
					# For rank 1 mtxiNdL[1][1], [1][4], [1][7],[1][10], [3][..] are wrong among more
					# These are only found wrong in loop C :: 132-134 >> externalToLocalMap which is 
					# Its wrong because indexing in the ReceiveList is all wrong. I used push! to add
					# elements to an existing set. Sets add these newly added values at the top instead
					# from at the bottom. Because of this indexing in receiveList, receiveEntryCount 
					# was all wrong and adding 162 in places it had to add 1.
					# SOLUTION 1: Make your own alternate set structure which maintains the chronological
					# 			  sense of what was added first to the set 
					# 			  type A : Instead of a set we make a dict with chronology stored as 
					# 			  		   as values to set vlaues(which will be keys)
					# 			  type B : Make 2 sets. Simulataneously add chronolgy and set values
					# 			  		   I am not so sure about this approach because I am not sure 
					# 			  		   if Julia will maintain these chronologies without a doubt  
					# For rank 0 all values match with C version 
					# mtxIndG is correct for both ranks for all positions

                    mtxIndL[i,j] = A.globalToLocalMap[curIndex]
###################################################
#		if A.geom.rank == 0 
#			open("j_same_rank_0.txt", "a") do f 
#				println(f,"$i, $j, $curIndex, $(A.globalToLocalMap[curIndex])") 
#			end
#		else
#			open("j_same_rank_1.txt", "a") do f 
#				println(f,"$i, $j, $curIndex, $(A.globalToLocalMap[curIndex])")
#			end
#		end
###################################################			
                else # If column index is not a row index, then it comes from another processor
                    mtxIndL[i,j] = externalToLocalMap[curIndex]
###################################################
#		if A.geom.rank == 0 
#			open("j_diff_rank_0.txt", "a") do f 
#				println(f,"$i, $j, $curIndex, $(externalToLocalMap[curIndex])") 
#			end
#		else
#			open("j_diff_rank_1.txt", "a") do f 
#				println(f,"$i, $j, $curIndex, $(externalToLocalMap[curIndex])")
#			end
#		end
###################################################			
                end
            end
       end

        # Store contents in our matrix struct
	if A.geom.rank==1
		open("rank1_elementsToSend", "a") do f
			DelimitedFiles.writedlm(f,elementsToSend)
		end 
	end
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

		if A.geom.rank == 0 
			open("j_mtxIndG_0.txt", "a") do f 
				println(f,"$(A.mtxIndG)") 
			end
		else
			open("j_mtxIndG_1.txt", "a") do f 
				println(f,"$(A.mtxIndG)") 
			end
		end
		if A.geom.rank == 0 
			open("j_mtxIndL_0.txt", "a") do f 
				println(f,"$(A.mtxIndL)") 
			end
		else
			open("j_mtxIndL_1.txt", "a") do f 
				println(f,"$(A.mtxIndL)") 
			end
		end

        @debug(" For rank $A.geom.rank of $A.geom.size number of neighbors $A.numberOfSendNeighbors")
        for i = 1:numberOfSendNeighbors
            @debug("     rank = ", A.geom.rank," neighbor = ",neighbors[i]," send/recv length = ", sendLength[i]/receiveLength[i],".")
            for j = 1:sendLength[i]
                @debug("       rank = ", A.geom.rank," elementsToSend[$j] =", elementsToSend[j],".")
            end
        end

    end # ! NO_MPI
    println("Out of SetupHalo_ref")
end
