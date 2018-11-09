
//@HEADER
// ***************************************************
//
// HPCG: High Performance Conjugate Gradient Benchmark
//
// Contact:
// Michael A. Heroux ( maherou@sandia.gov)
// Jack Dongarra     (dongarra@eecs.utk.edu)
// Piotr Luszczek    (luszczek@eecs.utk.edu)
//
// ***************************************************
//@HEADER

/*!
 @file SparseMatrix.hpp

 HPCG data structures for the sparse matrix
 */


include("Geometry.jl")
include("Vector.jl")
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
   mtxIndG # matrix indices as global values
   mtxIndL # matrix indices as local values
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
   Ac #Coarse grid matrix
   mgData #Pointer to the coarse level data for this fine matrix
   optimizationData  #pointer that can be used to store implementation-specific data

   numberOfSendNeighborsi::Int64 # number of neighboring processes that will be send local data
   totalToBeSent::Int64  # total number of entries to be sent
   elementsToSend #elements to send to neighboring processes
   neighbors #neighboring processes
   receiveLength # lenghts of messages received from neighboring processes
   sendLength # lenghts of messages sent to neighboring processes
   sendBuffer # send buffer for non-blocking sends
end

#=
  Initializes the known system matrix data structure members to 0.

  @param[in] A the known system matrix
=#
@inline function InitializeSparseMatrix(A, geom) 
  A.title = 0
  A.geom = geom
  A.totalNumberOfRows = 0
  A.totalNumberOfNonzeros = 0
  A.localNumberOfRows = 0
  A.localNumberOfColumns = 0
  A.localNumberOfNonzeros = 0
  A.nonzerosInRow = 0
  A.mtxIndG = 0
  A.mtxIndL = 0
  A.matrixValues = 0
  A.matrixDiagonal = 0

  # Optimization is ON by default. The code that switches it OFF is in the
  # functions that are meant to be optimized.
  A.isDotProductOptimized = true
  A.isSpmvOptimized       = true
  A.isMgOptimized      = true
  A.isWaxpbyOptimized     = true

  A.numberOfExternalValues = 0
  A.numberOfSendNeighbors = 0
  A.totalToBeSent = 0
  A.elementsToSend = 0
  A.neighbors = 0
  A.receiveLength = 0
  A.sendLength = 0
  A.sendBuffer = 0
  A.mgData = 0 #Fine-to-coarse grid transfer initially not defined.
  A.Ac =0
  return
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
	dv[i] = *(curDiagA[i])
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
    delete [] A.matrixValues[i]
    delete [] A.mtxIndG[i]
    delete [] A.mtxIndL[i]
  end

  delete [] A.matrixValues[1]
  delete [] A.mtxIndG[1]
  delete [] A.mtxIndL[1]

  if (A.title)
    delete [] A.title
  end
  if (A.nonzerosInRow)             
    delete [] A.nonzerosInRow
  end
  if (A.mtxIndG) 
    delete [] A.mtxIndG
  end
  if (A.mtxIndL) 
    delete [] A.mtxIndL
  end
  if (A.matrixValues) 
    delete [] A.matrixValues
  end
  if (A.matrixDiagonal)           
    delete [] A.matrixDiagonal
  end
  if (A.elementsToSend)       
    delete [] A.elementsToSend
  end
  if (A.neighbors)              
    delete [] A.neighbors
  end
  if (A.receiveLength)            
    delete [] A.receiveLength
  end
  if (A.sendLength)            
    delete [] A.sendLength
  end
  if (A.sendBuffer)            
    delete [] A.sendBuffer
  end
  if (A.geom!=0)  
    DeleteGeometry(*A.geom) 
    delete A.geom 
    A.geom = 0
  end
  if (A.Ac!=0) 
    DeleteMatrix(*A.Ac) 
    delete A.Ac A.Ac = 0 # Delete coarse matrix
  end
  if (A.mgData!=0) 
    DeleteMGData(*A.mgData) 
    delete A.mgData 
    A.mgData = 0} # Delete MG data
  end
  return
end

