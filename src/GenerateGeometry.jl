#=
  @file GenerateGeometry.cpp

  HPCG routine
=#

include("ComputeOptimalShapeXYZ.jl")

include("hpcg.jl")


#=
  Computes the factorization of the total number of processes into a
  3-dimensional process grid that is as close as possible to a cube. The
  quality of the factorization depends on the prime number structure of the
  total number of processes. It then stores this decompostion together with the
  parallel parameters of the run in the geometry data structure.

  @param[in]  size total number of MPI processes
  @param[in]  rank this process' rank among other MPI processes
  @param[in]  numThreads number of OpenMP threads in this process
  @param[in]  pz z-dimension processor ID where second zone of nz values start
  @param[in]  nx, ny, nz number of grid points for each local block in the x, y, and z dimensions, respectively
=#
function generate_geometry!(size, rank, numThreads,pz, zl, zu,
  			nx, ny, nz, npx, npy, npz)
  if npx * npy * npz <= 0 || npx * npy * npz > size
    if MPI.Initialized == false
	size = 1
    end
    npx, npy, npz = compute_optimal_shape_xyz(size, npx, npy, npz)
  end

  @debug("npx = $npx, npy = $npy, npz = $npz, comparing factor = $(npx*npy*npz)")
  @debug("npz = $npz")

  npartz    = 0

  if pz == 0 # No variation in nz sizes
    npartz       = 1
    partz_ids    = Array{Int64,1}(undef,1) # WHAT
    partz_nz     = Array{Int64,1}(undef,1) #WHAT
    partz_ids[1] = npz  #CHECK INDICES HERE
    partz_nz[1]  = nz
  else 
    npartz       = 2
    partz_ids    = Array{Int64}(undef,2)
    partz_ids[1] = pz
    partz_ids[2] = npz
    partz_nz     = Array{Int64}(undef,2)
    partz_nz[1]  = zl
    partz_nz[2]  = zu
  end

#  partz_ids[npartz-1] = npz # The last element of this array is always npz
  ipartz_ids = 0
  
  for i=1 :npartz 
#    @assert(ipartz_ids<partz_ids[i])  # Make sure that z partitioning is consistent with computed npz value
    ipartz_ids = partz_ids[i]
  end

  # Now compute this process's indices in the 3D cube
  ipz = rank ÷ (npx*npy)
  ipy = (rank-ipz*npx*npy)÷npx
  ipx = rank%npx # will gice division error because npx is zero 

  @debug("Generate Geometry: npx, npy, npz, nx, ny, nz, ipx, ipy, ipz => $npx, $npy, $npz, $nx, $ny, $nz, $ipx, $ipy, $ipz")

  if rank == 0
    @debug("size = $size\n
        nx  = $nx\n
        ny  = $ny\n
        nz  = $nz\n
        npx = $npx\n
        npy = $npy\n
        npz = $npz\n")

     @debug("For rank = $rank\n
      ipx = $ipx\n
      ipy = $ipy\n
      ipz = $ipz\n")

  end

  @assert(size>=npx*npy*npz)

# These values should be defined to take into account changes in nx, ny, nz values
# due to variable local grid sizes
  gnx = npx*nx
  gny = npy*ny
  gnz = 0

  # We now permit varying values for nz for any nx-by-ny plane of MPI processes.
  # npartz is the number of different groups of nx-by-ny groups of processes.
  # partz_ids is an array of length npartz where each value indicates the z process of the last process in the ith nx-by-ny group.
  # partz_nz is an array of length npartz containing the value of nz for the ith group.

  #        With no variation, npartz = 1, partz_ids[0] = npz, partz_nz[0] = nz

  ipartz_ids = 0

  for i=1:npartz
    ipartz_ids = partz_ids[i] - ipartz_ids
    gnz += partz_nz[i]*ipartz_ids
  end

  giz0       = 0
  ipartz_ids = 0

  for i=1:npartz 
    ipart_nz = partz_nz[i]
    if (ipz < partz_ids[i]) 
      giz0 += (ipz-ipartz_ids)*ipart_nz
      break
     else 
      ipartz_ids = partz_ids[i]
      giz0 += ipartz_ids*ipart_nz
    end
  end

  gix0 = ipx*nx
  giy0 = ipy*ny


  geom = Geometry(size, rank, numThreads, 
                  nx, ny, nz, npx, npy, npz, 
                  pz, npartz, partz_ids, partz_nz, 
                  ipx, ipy, ipz, 
                  gnx, gny, gnz, 
                  gix0, giy0, giz0)

  return geom

end
