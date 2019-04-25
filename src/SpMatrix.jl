#=
 @file SparseMatrix.hpp

 HPCG data structures for the sparse matrix
=#


include("Geometry.jl")
include("MGData.jl")

mutable struct HPCGSparseMatrix

    is_dot_prod_optimized::Bool
    is_spmv_optimized::Bool
    is_mg_optimized::Bool
    is_waxpby_optimized::Bool

    geom::Geometry                            # geometry associated with this matrix

    title::String                             # name of the sparse matrix
    totalNumberOfRows::Int64                  # total number of matrix rows across all processes
    totalNumberOfNonzeros::Int64              # total number of matrix nonzeros across all processes
    localNumberOfRows::Int64                  # number of rows local to this process
    localNumberOfColumns::Int64               # number of columns local to this process
    localNumberOfNonzeros::Int64              # number of nonzeros local to this process
    nonzerosInRow::Array{Any}                             # The number of nonzeros in a row will always be 27 or fewer
    mtxIndG ::Array{Array{Int64,1}}           # matrix indices as global values
    mtxIndL ::Array{Array{Int64,1}}           # matrix indices as local value
    matrixValues :: Array{Array{Float64,1}}   # values of matrix entries
    matrixDiagonal :: Array{Array{Float64,1}} # values of matrix diagonal entries
    localToGlobalMap::Dict    # local-to-global mapping
    globalToLocalMap::Dict    # global-to-local mapping


    numberOfExternalValues::Int64
    numberOfSendNeighbors::Int64 # number of neighboring processes that will be send local data
    totalToBeSent::Int64         # total number of entries to be sent
    elementsToSend::Array{Int64} # elements to send to neighboring processes
    neighbors::Array{Int64}      # neighboring processes
    receiveLength::Array{Int64}  # lenghts of messages received from neighboring processes
    sendLength::Array{Int64}     # lenghts of messages sent to neighboring processes
    sendBuffer::Array{Float64}    # send buffer for non-blocking sends

    #=
    This is for storing optimized data structres created in OptimizeProblem and
    used inside optimized compute_spmv().
    =#

    Ac::HPCGSparseMatrix  # Coarse grid matrix
    mgData::MGData       # Pointer to the coarse level data for this fine matrix


    function HPCGSparseMatrix(dpopt, spmvopt, mgopt, waxpbyopt, g,
                              ttl, tnrows, tnnz, lnrows, lncols, lnnz,
                              nzinrow, mtxindg, mtxindl, matvals, matdiag, 
                              l2gmap, g2lmap,
                              nextvals, nsendneighbor, totaltosend,
                              elmtosend, neighb, rcvlen, sendlen, sendbuf,
                              mgd)

        x = new()

        x.is_dot_prod_optimized  = dpopt
        x.is_spmv_optimized      = spmvopt
        x.is_mg_optimized        = mgopt
        x.is_waxpby_optimized    = waxpbyopt

        x.geom                   = g

        x.title                  = ttl
        x.totalNumberOfRows      = tnrows
        x.totalNumberOfNonzeros  = tnnz
        x.localNumberOfRows      = lnrows
        x.localNumberOfColumns   = lncols
        x.localNumberOfNonzeros  = lnnz
        x.nonzerosInRow          = nzinrow
        x.mtxIndG                = mtxindg
        x.mtxIndL                = mtxindl
        x.matrixValues           = matvals
        x.matrixDiagonal         = matdiag
        x.localToGlobalMap       = l2gmap
        x.globalToLocalMap       = g2lmap

        x.numberOfExternalValues = nextvals
        x.numberOfSendNeighbors  = nsendneighbor
        x.totalToBeSent          = totaltosend
        x.elementsToSend         = elmtosend
        x.neighbors              = neighb
        x.receiveLength          = rcvlen
        x.sendLength             = sendlen
        x.sendBuffer             = sendbuf

        x.Ac                     = new()
        x.mgData                 = mgd

        return x

    end 

end


function HPCGSparseMatrix(dprod_opt, spmb_opt, mg_opt, waxpby_opt, geom)
    return HPCGSparseMatrix(dprod_opt, spmb_opt, mg_opt, waxpby_opt, geom,
                            "", 0, 0, 0, 0, 0, [], [], [], [], [], Dict(), Dict(),
                            0, 0, 0, [], [], [], [], [], MGData())

end

#=
  Initializes the known system matrix data structure members to 0.

  @param[in] A the known system matrix
=#
function initialize_sparse_matrix(geom::Geometry) 
  A = HPCGSparseMatrix(true, true , true, true, geom)
  return A
end

#=
  Copy values from matrix diagonal into user-provided vector.

  @param[in] A the known system matrix.
  @param[inout] diagonal  Vector of diagonal values (must be allocated before call to this function).
=#
@inline function copy_matrix_diagonal(A) 
    cur_diag_a = A.matrixDiagonal
    dv         = Vector{Int64}(undef, A.localNumberOfRows)
    @assert(A.localNumberOfRows == length(dv))
    dv         = cur_diag_a
    return dv
end

#=
  Replace specified matrix diagonal value.

  @param[inout] A The system matrix.
  @param[in] diagonal  Vector of diagonal values that will replace existing matrix diagonal values.
=#
@inline function replace_matrix_diagonal!(A, diag) 
    @assert(A.localNumberOfRows==length(diag))
    A.matrixDiagonal = diag
end
