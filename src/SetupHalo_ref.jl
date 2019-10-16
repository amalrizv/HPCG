#=
 @file SetupHalo_ref.jl

 HPCG routine
=#
using DataStructures
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

		sendList           = Dict{Int64, OrderedSet{Int64}}()
		receiveList        = Dict{Int64, OrderedSet{Int64}}()
        externalToLocalMap = Dict{Int64, Int64}()
		#BUG_INFO #RZV
		# The main problem is with the Set, OrderedSet data Structure that is used here
		# TRIAL 1 : I tried using a Set first but since it did not retain order we switched to
		# TRIAL 2 : a custom Dict which was similar to the PriorityQueue structure in DataStructures.
		# TRIAL 3 : To remove unnneccessary complications (TRIAL_2) we just arranged 
		# 			the elements of the Set in an ascending order using sort(collect...) 
		# 			but that was a wrong approach since elements of Set are not in 
		# 			ascending order but just FIFO.
		# TRIAL_4 : So I switched to using an OrderedSet.
		# 			An OrderedSet is a LIFO Set. Operated by using push! and pop!
		# 			when the values are collected using collect(some_ordered_set)
		# 			the FIFO order is maintained. [Output evidence  in set_output_1]
		#			+
		#			when the values are iterated through the Ordered Set then also 
		#			a FIFO order is maintained.
		#			= Even after maintaining FIFO nature in iteration and collect (in REPL and ord.jl) 
		#			  it is seen that when the code is assigning externalToLocalMap
		#			  for keys that are either being iterated via ordered_set_values or 
		#			  collected_ordered_set_values or baffilingly even 
		#			  reverse_collected_ordered_set_values. Some stubborn unknown condition
		#			  in the code is making these values being read in LIFO order.
		#
		#
		#
		for ranks = 0: A.geom.size-1
			receiveList[ranks] = OrderedSet{Int64}()
			sendList[ranks]     = OrderedSet{Int64}()
		end

        #  TODO: With proper critical and atomic regions, this loop could be threaded, but not attempting it at this time
		#
        for i = 1:localNumberOfRows 
            currentGlobalRow = A.localToGlobalMap[i] 
            for j = 1:nonzerosInRow[i]
				curIndex            = mtxIndG[i,j]
                rankIdOfColumnEntry = compute_rank_of_matrix_row(A.geom, curIndex)

###################################################
#CORRECT
				if A.geom.rank == 0 
					open("rankIdOfColumnEntry_0.txt", "a") do f 
						println(f, "rankIdOfColumnEntry is $rankIdOfColumnEntry for curIndex(mtxIndG)[$i][$j] value $curIndex")
					end
				else
					open("rankIdOfColumnEntry_1.txt", "a") do f 
						println(f, "rankIdOfColumnEntry is $rankIdOfColumnEntry for curIndex(mtxIndG)[$i][$j] value $curIndex")
					end
#CORRECT
				end
###################################################

            	if  A.geom.rank != rankIdOfColumnEntry # If column index is not a row index, then it comes from another processor
					
		    		push!(receiveList[rankIdOfColumnEntry], curIndex)

                	push!(sendList[rankIdOfColumnEntry], currentGlobalRow) 
					
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
        # @debug("totalToBeSent = $totalToBeSent totalToBeReceived = $totalToBeReceived")
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
        neighborCount     = 0
		receiveEntryCount = 0 
		sendEntryCount	  = 0


        for (k,v) in receiveList 
            neighborId                   = k # rank of current neighbor we are processing
            neighbors[neighborCount+1]     = neighborId # store rank ID of current neighbor
            receiveLength[neighborCount+1] = length(v)
            sendLength[neighborCount+1]    = length(sendList[neighborId]) # Get count if sends/receives
###################################################
#=
			n_rcv_id = reverse!(collect(v)) 							# n_rcv_id is an array of  set elements
																# no element repeated

			n_snd_id = reverse!(collect(sendList[neighborId]))		# n_snd_id is an array of set elements
																# no element repeated
=#					
			if A.geom.rank == 0
				open("set_output_0.txt", "a") do f 
					println(f, "receiveList[k] \n $(receiveList[k])")
					println(f,"collect(receiveList[k] \n $(collect(receiveList[k]))")
					println(f, "reverse!(collect(receiveList[k] \n $(reverse!(collect(receiveList[k])))")
				end
			else
				open("set_output_1.txt", "a") do f 
					println(f, "receiveList[k] \n $(receiveList[k])")
					println(f,"collect(receiveList[k] \n $(collect(receiveList[k]))")
					println(f, "reverse!(collect(receiveList[k] \n $(reverse!(collect(receiveList[k])))")
				end
			end
			for x in sort(collect(receiveList[k]))
		

###################################################
				if A.geom.rank == 0 
					open("j_e2lmapping_0.txt", "a") do f 
						println(f,"externalToLocalMap[$x] = $localNumberOfRows + $receiveEntryCount +1") 
					end
				else
					open("j_e2lmapping_1.txt", "a") do f 
						println(f,"externalToLocalMap[$x] = $localNumberOfRows + $receiveEntryCount +1") 
					end

				end

				
###################################################
				#RZV #BUG_INFO
				# test fucntion file  = ord.jl
				# Dict(externalToLocalMap) does not maintain the order in which elements are sent
				# to receiveList[rank_id] but that is not to say that the mapping would be wrong.
				#
				# My assumption was that tha order might change but the mapping will be correct
				#  
				# Is this behaviour seen in both processes ?  Yes
				#
				# Is this behaviour tested in Julia REPL 
				# 			or a seperate testing function ?  Yes the order of dict changes but the mapping is preserved
				#
				# So behaviour of Test function and this code are not behaving the same way.

		# BUG_INFO #RZV
		# Because of all this outputs recorded for $externalToLocalMap[curIndex] for all ranks is wrong
		# Because of all this outputs recorded for j_e2lmapping_0&1 j_diff_rank_0&1 and j_mtxIndL_0&1 for all ranks is wrong
				externalToLocalMap[x] = localNumberOfRows + receiveEntryCount + 1 # The remote columns are indexed at end of internals
				receiveEntryCount    += 1
            end

			for x in sort(collect(sendList[k]))

###################################################
				if A.geom.rank == 0 
					open("els_2_send_rank_0.txt", "a") do f 
						println(f,"elementsToSend[$sendEntryCount +1] = $(A.globalToLocalMap[x])") 
					end
				else
					open("els_2_send_rank_1.txt", "a") do f 
						println(f,"elementsToSend[$sendEntryCount +1] = $(A.globalToLocalMap[x])") 
					end
				end
###################################################

			elementsToSend[sendEntryCount+1] = A.globalToLocalMap[x] # store local ids of entry to send
				sendEntryCount   += 1
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

###################################################
		if A.geom.rank == 0 
			open("j_same_rank_0.txt", "a") do f 
				println(f,"$i, $j, $curIndex, $(A.globalToLocalMap[curIndex])") 
			end
		else
			open("j_same_rank_1.txt", "a") do f 
				println(f,"$i, $j, $curIndex, $(A.globalToLocalMap[curIndex])")
			end
		end
###################################################			
                else # If column index is not a row index, then it comes from another processor

		# BUG_INFO #RZV
		# Because of all this outputs recorded for $externalToLocalMap[curIndex] for all ranks is wrong
		# Because of all this outputs recorded for j_e2lmapping_0&1 j_diff_rank_0&1 and j_mtxIndL_0&1 for all ranks is wrong

		
                    mtxIndL[i,j] = externalToLocalMap[curIndex]
###################################################
		if A.geom.rank == 0 
			open("j_diff_rank_0.txt", "a") do f 
				println(f,"$i, $j, $curIndex, $(externalToLocalMap[curIndex])") 
			end
		else
			open("j_diff_rank_1.txt", "a") do f 
				println(f,"$i, $j, $curIndex, $(externalToLocalMap[curIndex])")
			end

		end
###################################################			
#
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
		# BUG_INFO #RZV
		# Because of all this outputs recorded for mtxIndL for all ranks is wrong
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
