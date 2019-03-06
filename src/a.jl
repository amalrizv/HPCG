using MPI
include("hpcg.jl")
function GenerateProblem_ref(A, b,x, xexact) 

  nx = A.geom.nx
  ny = A.geom.ny
  nz = A.geom.nz
  gnx = A.geom.gnx
  gny = A.geom.gny
  gnz = A.geom.gnz
  gix0 = A.geom.gix0
  giy0 = A.geom.giy0
  giz0 = A.geom.giz0

  localNumberOfRows = nx*ny*nz # This is the size of our subblock
  @assert(localNumberOfRows>0) # Throw an exception of the number of rows is less than zero (can happen if int overflow)
  numberOfNonzerosPerRow = 27 # We are approximating a 27-point finite element/volume/difference 3D stencil

  totalNumberOfRows = gnx*gny*gnz # Total number of grid points in mesh
  @assert(totalNumberOfRows>0) # Throw an exception of the number of rows is less than zero (can happen if int overflow)


  nonzerosInRow = Array{Any}(undef,localNumberOfRows)
  mtxIndG = Array{Int64}(undef, localNumberOfRows)
  mtxIndL =Array{Int64}(undef,localNumberOfRows)
  matrixValues = Array{Float64}(undef,localNumberOfRows)
  matrixDiagonal = Array{Float64}(undef,localNumberOfRows)

  if b!=0 
	b  = Vector(localNumberOfRows)
  end
  if x!=0
	x =  Vector(localNumberOfRows)
  end
  if xexact!=0
	xexact = Vector(localNumberOfRows)
  end
  bv = zeros(localNumberOfRows)
  xv = zeros(localNumberOfRows)
  xexactv = zeros(localNumberOfRows)
  if b!=0
	 bv = b # Only compute exact solution if requested
  end
  if x!=0
	 xv = x # Only compute exact solution if requested
  end
  if xexact!=0
	 xexactv = xexact # Only compute exact solution if requested
  end

  resize!(A.localToGlobalMap, localNumberOfRows)


  for i=1:localNumberOfRows 
    	matrixValues[i] = 0
    	matrixDiagonal[i] = 0
    	mtxIndG[i] = 0
    	mtxIndL[i] = 0
  end
  for i=1:localNumberOfRows 
    	mtxIndL[i] = Array{Float64}(undef,numberOfNonzerosPerRow)
  end
  for i=1:localNumberOfRows 
    	matrixValues[i] = Array{Float64}(undef,numberOfNonzerosPerRow)
  end
  for i=1:localNumberOfRows
   	mtxIndG[i] = Array{Float64}(undef,numberOfNonzerosPerRow)
  end
  mtxIndL[1] = Array{Float64}(undef,localNumberOfRows * numberOfNonzerosPerRow)
  matrixValues[1] = Array{Float64}(undef,localNumberOfRows * numberOfNonzerosPerRow)
  mtxIndG[1] = Array{Float64}(undef,localNumberOfRows * numberOfNonzerosPerRow)

   localNumberOfNonzeros = 0
  return
end
