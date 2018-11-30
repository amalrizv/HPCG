include("GenerateCoarseProblem.jl")
include("GenerateGeometry.jl")
include("GenerateProblem.jl")
include("SetupHalo.jl")

#=
  Routine to construct a prolongation/restriction operator for a given fine grid matrix
  solution (as computed by a direct solver).

  @param[inout]  Af - The known system matrix, on output its coarse operator, fine-to-coarse operator and auxiliary vectors will be defined.

  Note that the matrix Af is considered const because the attributes we are modifying are declared as mutable.
=#

function GenerateCoarseProblem(const Af) 

  # Make local copies of geometry information.  Use global_int_t since the RHS products in the calculations
  # below may result in global range values.
   nxf = Af.geom.nx
   nyf = Af.geom.ny
   nzf = Af.geom.nz

  nxc = Int64
  nyc = Int64
  nzc = Int64 #Coarse nx, ny, nz
  @assert(nxf%2==0) 
  @assert(nyf%2==0) 
  @assert(nzf%2==0) # Need fine grid dimensions to be divisible by 2
  nxc = nxf/2 nyc = nyf/2 nzc = nzf/2
  f2cOperator = Array{Int64}(undef,Af.localNumberOfRows)
   localNumberOfRows = nxc*nyc*nzc # This is the size of our subblock
  # If this @assert fails, it most likely means that the local_int_t is set to int and should be set to long long
  @assert(localNumberOfRows>0) # Throw an exception of the number of rows is less than zero (can happen if "int" overflows)

  # Use a parallel loop to do initial assignment:
  # distributes the physical placement of arrays of pointers across the memory system
  for  i=1:localNumberOfRows 
    f2cOperator[i] = 0
  end


  # TODO:  This triply nested loop could be flattened or use nested parallelism
  for  izc=0:nzc 
     izf = 2*izc
    for iyc=0:nyc
       iyf = 2*iyc
      for ixc=0:nxc 
         ixf = 2*ixc
         currentCoarseRow = izc*nxc*nyc+iyc*nxc+ixc
         currentFineRow = izf*nxf*nyf+iyf*nxf+ixf
        f2cOperator[currentCoarseRow] = currentFineRow
      end # end iy loop
    end # end even iz if statement
  end # end iz loop

  # Construct the geometry and linear system
  geomc = Geometry
  zlc = 0 # Coarsen nz for the lower block in the z processor dimension
  zuc = 0 # Coarsen nz for the upper block in the z processor dimension
  pz = Af.geom.pz
  if pz>0
    zlc = Af.geom.partz_nz[0]/2 # Coarsen nz for the lower block in the z processor dimension
    zuc = Af.geom.partz_nz[1]/2 # Coarsen nz for the upper block in the z processor dimension
  end
  GenerateGeometry(Af.geom.size, Af.geom.rank, Af.geom.numThreads, Af.geom.pz, zlc, zuc, nxc, nyc, nzc, Af.geom.npx, Af.geom.npy, Af.geom.npz, geomc)

  Ac = SpMatrix
  InitializeSparseMatrix(Ac, geomc)
  GenerateProblem(Ac, 0, 0, 0)
  SetupHalo(Ac)
  rc = new Vector
  xc = new Vector
  Axf = new Vector
  InitializeVector(rc, Ac.localNumberOfRows)
  InitializeVector(xc, Ac.localNumberOfColumns)
  InitializeVector(Axf, Af.localNumberOfColumns)
  Af.Ac = Ac
  mgData = MGData
  InitializeMGData(f2cOperator, rc, xc, Axf, mgData)
  Af.mgData = mgData

  return
end
