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
  #Need fine grid dimensions to be divisible by 2
  @assert(nxf%2==0) 
  @assert(nyf%2==0) 
  @assert(nzf%2==0) 
  nxc = cld(nxf,2) 
  nyc = cld(nyf,2)
  nzc = cld(nzf,2)
  f2cOperator = Array{Int64}(undef,A.localNumberOfRows)
   localNumberOfRows = nxc*nyc*nzc # This is the size of our subblock
  # If this @assert fails, it most likely means that the local_int_t is set to int and should be set to long long
  @assert(localNumberOfRows>0) # Throw an exception of the number of rows is less than zero (can happen if "int" overflows)

  # Use a parallel loop to do initial assignment:")
  # distributes the physical placement of arrays of pointers across the memory system
  for  i=1:localNumberOfRows 
    f2cOperator[i] = 0
  end


  # TODO:  This triply nested loop could be flattened or use nested parallelism
  for  izc=1:nzc 
     izf = 2*(izc-1)+1
    for iyc=1:nyc
       iyf = 2*(iyc-1)+1
      for ixc=1:nxc 
         ixf = 2*(ixc-1)+1
         currentCoarseRow = izc*nxc*nyc+iyc*nxc+ixc
         currentFineRow = izf*nxf*nyf+iyf*nxf+ixf
        f2cOperator[currentCoarseRow] = currentFineRow
      end # end iy loop
    end # end even iz if statement
  end # end iz loop

  # Construct the geometry and linear system")
  zlc = 0 # Coarsen nz for the lower block in the z processor dimension
  zuc = 0 # Coarsen nz for the upper block in the z processor dimension
  pz  = A.geom.pz

  if pz>0
    zlc = cld(A.geom.partz_nz[1],2) # Coarsen nz for the lower block in the z processor dimension
    zuc = cld(A.geom.partz_nz[2],2) # Coarsen nz for the upper block in the z processor dimension
  end

  geomc = generate_geometry(A.geom.size, A.geom.rank, A.geom.numThreads, A.geom.pz, zlc, zuc, nxc, nyc, nzc, A.geom.npx, A.geom.npy, A.geom.npz)

  Ac               = initialize_sparse_matrix(geomc) 	
  ret1, ret2, ret3 =  generate_problem!(Ac)		
  setup_halo!(Ac)
  rc          = Vector{Int64}(undef, localNumberOfRows)
  xc          = Vector{Int64}(undef, localNumberOfRows)
  Axf         = Vector{Int64}(undef, localNumberOfRows)
  mgd::MGData = InitializeMGData(f2cOperator, rc, xc, Axf)
  A.Ac        = Ac 
  A.MGData    = mgd	#sp_coarse structure where Ac is sp_anx structure
end

