#=
 @file SparseMatrix.hpp

 HPCG data structures for the sparse matrix
=#


include("Geometry.jl")
include("MGData.jl")

mutable struct SpMatrix
  title::String #name of the sparse matrix
  geom::Geometry #geometry associated with this matrix
   totalNumberOfRows::Int64 #total number of matrix rows across all processes
   totalNumberOfNonzeros::Int64  #total number of matrix nonzeros across all processes
   localNumberOfRows::Int64  # number of rows local to this process
   localNumberOfColumns::Int64   #number of columns local to this process
   localNumberOfNonzeros::Int64   # number of nonzeros local to this process
   nonzerosInRow  # The number of nonzeros in a row will always be 27 or fewer
   mtxIndG ::Array{Int64,2}# matrix indices as global values
   mtxIndL ::Array{Int64,2}# matrix indices as local values
   matrixValues # values of matrix entries
   matrixDiagonal # values of matrix diagonal entries
   globalToLocalMap #global-to-local mapping
   localToGlobalMap #local-to-global mapping
   isDotProductOptimized
   isSpmvOptimized
   isMgOptimized
   isWaxpbyOptimized
  #=
   This is for storing optimized data structres created in OptimizeProblem and
   used inside optimized ComputeSPMV().
  =#

   numberOfSendNeighbors::Int64 # number of neighboring processes that will be send local data
   totalToBeSent::Int64  # total number of entries to be sent
   elementsToSend #elements to send to neighboring processes
   neighbors #neighboring processes
   receiveLength # lenghts of messages received from neighboring processes
   sendLength # lenghts of messages sent to neighboring processes
   sendBuffer # send buffer for non-blocking sends
   mgData #Pointer to the coarse level data for this fine matrix
   Ac #Coarse grid matrix

end

#=
  Initializes the known system matrix data structure members to 0.

  @param[in] A the known system matrix
=#
function InitializeSparseMatrix(A::Type{SpMatrix} , geom::Geometry) 
   #Optimization data not inlcuded in struct
  A = SpMatrix("", geom,0,0,0,0,0,0,0,0,0,0,0,0, true, true, true,true, 0,0,0,0,0,0,0,0,0)
  return A
end

#=
  Copy values from matrix diagonal into user-provided vector.

  @param[in] A the known system matrix.
  @param[inout] diagonal  Vector of diagonal values (must be allocated before call to this function).
=#
@inline function CopyMatrixDiagonal(A, diagonal) 
    curDiagA = A.matrixDiagonal
    dv = diagonal.values
    @assert(A.localNumberOfRows==diagonal.localLength)
    for i=1:A.localNumberOfRows 
	dv[i] = curDiagA[i]
    end
  return
end
#=
  Replace specified matrix diagonal value.

  @param[inout] A The system matrix.
  @param[in] diagonal  Vector of diagonal values that will replace existing matrix diagonal values.
=#
@inline function ReplaceMatrixDiagonal(A, diagonal) 
    curDiagA = A.matrixDiagonal
    dv = diagonal.values
    @assert(A.localNumberOfRows==diagonal.localLength)
    for i=1:A.localNumberOfRows
	 curDiagA[i] = dv[i]
    end
  return
end
#=
  Deallocates the members of the data structure of the known system matrix provided they are not 0.

  @param[in] A the known system matrix
=#
@inline function  DeleteMatrix(A) 

  for i = 1:A.localNumberOfRows 
    A.matrixValues[i] = nothing
    A.mtxIndG[i] = nothing
    A.mtxIndL[i] = nothing
  end

  A.matrixValues[1] = nothing
  A.mtxIndG[1] = nothing
  A.mtxIndL[1] = nothing

  if (A.title)
    A.titl = nothing
  end
  if (A.nonzerosInRow)             
    A.nonzerosInRow = nothing
  end
  if (A.mtxIndG) 
    A.mtxIndG = nothing
  end
  if (A.mtxIndL) 
    A.mtxIndL = nothing
  end
  if (A.matrixValues) 
    A.matrixValues = nothing
  end
  if (A.matrixDiagonal)           
    A.matrixDiagonal = nothing
  end
  if (A.elementsToSend)       
    A.elementsToSend = nothing
  end
  if (A.neighbors)              
    A.neighbors = nothing
  end
  if (A.receiveLength)            
    A.receiveLength = nothing
  end
  if (A.sendLength)            
    A.sendLength = nothing
  end
  if (A.sendBuffer)            
    A.sendBuffer = nothing
  end
  if (A.geom!=0)  
    DeleteGeometry(A.geom) 
    A.geom  = nothing
  end
  if (A.Ac!=0) 
    DeleteMatrix(A.Ac) 
    A.Ac = nothing # Delete coarse matrix
  end
  if (A.mgData!=0) 
    DeleteMGData(A.mgData) 
    A.mgData = nothing 
    # Delete MG data
  end
  return
end

