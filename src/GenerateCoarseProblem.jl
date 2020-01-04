include("GenerateGeometry.jl")
include("MGData.jl")
include("GenerateProblem.jl")
include("SetupHalo.jl")

#=
  Routine to construct a prolongation/restriction operator for a given fine grid matrix
  solution (as computed by a direct solver).

  @param[inout]  A - The known system matrix, on output its coarse operator, fine-to-coarse operator and auxiliary vectors will be defined.

  Note that the matrix A is considered const because the attributes we are modifying are declared as mutable.
=#

function generate_coarse_problem!(A)

    # Make local copies of geometry information.  Use global_int_t since the RHS products in the calculations
    # below may result in global range values.
    nxf = A.geom.nx
    nyf = A.geom.ny
    nzf = A.geom.nz

    nxc = Int64
    nyc = Int64
    nzc = Int64 #Coarse nx, ny, nz

    # Need fine grid dimensions to be divisible by 2
    @assert(nxf%2==0) 
    @assert(nyf%2==0) 
    @assert(nzf%2==0) 

    nxc = div(nxf, 2) 
    nyc = div(nyf, 2)
    nzc = div(nzf, 2)

    f2c_operator = Array{Int64}(undef,A.localNumberOfRows)
    local_num_rows = nxc*nyc*nzc # This is the size of our subblock:
    # If this @assert fails, it most likely means that the local_int_t is set to int and should be set to long long
    @assert(local_num_rows>0) # Throw an exception of the number of rows is less than zero (can happen if "int" overflows)

    # NOTE: originally omp_parallel_for
    # Use a parallel loop to do initial assignment:")
    # distributes the physical placement of arrays of pointers across the memory system
    fill!(f2c_operator,0)
#    for  i=1:local_num_rows 
#        f2c_operator[i] = 0
#    end


    # TODO:  This triply nested loop could be flattened or use nested parallelism
    #
    # KCH TODO: the indices below are probably incorrect
    for  izc=1:nzc 
        izf = 2*(izc-1)+1
        for iyc=1:nyc
            iyf = 2*(iyc-1)+1
            for ixc=1:nxc 
                ixf = 2*(ixc-1)+1
                currentCoarseRow = (izc-1)*nxc*nyc+(iyc-1)*nxc+(ixc-1)+1
                currentFineRow = (izf-1)*nxf*nyf+(iyf-1)*nxf+(ixf-1)+1
                f2c_operator[currentCoarseRow] = currentFineRow
            end # end iy loop
        end # end even iz if statement
    end # end iz loop

    # Construct the geometry and linear system")
    zlc = 0 # Coarsen nz for the lower block in the z processor dimension
    zuc = 0 # Coarsen nz for the upper block in the z processor dimension
    pz  = A.geom.pz

    if pz>0
        zlc = div(A.geom.partz_nz[1],2) # Coarsen nz for the lower block in the z processor dimension
        zuc = div(A.geom.partz_nz[2],2) # Coarsen nz for the upper block in the z processor dimension
    end

    geomc = generate_geometry!(A.geom.size, A.geom.rank, A.geom.numThreads, A.geom.pz, zlc, zuc, nxc, nyc, nzc, A.geom.npx, A.geom.npy, A.geom.npz)

    Ac               = initialize_sparse_matrix(geomc)
    ret2, ret3, ret4 = generate_problem!(Ac)

    setup_halo!(Ac)

    rc          = zeros(Ac.localNumberOfRows)
    xc          = zeros(Ac.localNumberOfColumns)
    
	Axf         = zeros(A.localNumberOfColumns)
    #rc          = Vector{Float64}(undef, Ac.localNumberOfRows)
    #xc          = Vector{Float64}(undef, Ac.localNumberOfColumns)
    #Axf         = Vector{Float64}(undef, A.localNumberOfColumns)
    
    mgd::MGData = init_mg_data(f2c_operator, rc, xc, Axf)

    A.Ac        = Ac
    A.mgData    = mgd
end

