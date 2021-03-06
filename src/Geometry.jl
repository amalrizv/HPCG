#=
 @file Geometry.jl

 HPCG data structure for problem geometry
=#

#= 
  This is a data structure to contain all processor geometry information
=#
struct Geometry
  size::Int64       # Number of MPI processes
  rank::Int64       # This process' rank in the range [0 to size - 1]
  numThreads::Int64 # This process' number of threads
  nx::Int64         # Number of x-direction grid points for each local subdomain
  ny::Int64         # !< Number of y-direction grid points for each local subdomain
  nz::Int64         # !< Number of z-direction grid points for each local subdomain
  npx::Int64        # !< Number of processors in x-direction
  npy::Int64        # !< Number of processors in y-direction
  npz::Int64        # !< Number of processors in z-direction
  pz::Int64         # !< partition ID of z-dimension process that starts the second region of nz values
  npartz::Int64     # !< Number of partitions with varying nz values
  partz_ids         # !< Array of partition ids of processor in z-direction where new value of nz starts (valid values are 1 to npz)
  partz_nz          # !< Array of length npartz containing the nz values for each partition
  ipx::Int64        # !< Current rank's x location in the npx by npy by npz processor grid
  ipy::Int64        # !< Current rank's y location in the npx by npy by npz processor grid
  ipz::Int64        # !< Current rank's z location in the npx by npy by npz processor grid
  gnx::Int64        # !< Global number of x-direction grid points
  gny::Int64        # !< Global number of y-direction grid points
  gnz::Int64        # !< Global number of z-direction grid points
  gix0::Int64       # !< Base global x index for this rank in the npx by npy by npz processor grid
  giy0::Int64       # !< Base global y index for this rank in the npx by npy by npz processor grid
  giz0::Int64       # !< Base global z index for this rank in the npx by npy by npz processor grid
end

#=
  Returns the rank of the MPI process that is assigned the global row index
  given as the input argument.

  @param[in] geom  The description of the problem's geometry.
  @param[in] index The global row index

  @return Returns the MPI rank of the process assigned the row
=#
function compute_rank_of_matrix_row(geom::Geometry, index::Int64) 
    gnx = geom.gnx
    gny = geom.gny
    iz = (index-1)÷(gny*gnx)
    iy = ((index-1)-iz*gny*gnx)÷gnx
    ix = (index-1)%gnx
    # We now permit varying values for nz for any nx-by-ny plane of MPI processes.
    # npartz is the number of different groups of nx-by-ny groups of processes.
    # partz_ids is an array of length npartz where each value indicates the z process of the last process in the ith nx-by-ny group.
    # partz_nz is an array of length npartz containing the value of nz for the ith group.

    #        With no variation, npartz = 1, partz_ids[0] = npz, partz_nz[0] = nz

    ipz        = 0
    ipartz_ids = 0

    for i=1:geom.npartz 
        ipart_nz = geom.partz_nz[i]
        ipartz_ids = geom.partz_ids[i] - ipartz_ids
        if iz<= ipart_nz*ipartz_ids
            ipz += iz÷ipart_nz
            break
        else 
            ipz += ipartz_ids
            iz -= ipart_nz*ipartz_ids
        end #if loop

    end #for loop

    #  ipz = iz/geom.nz
    ipy  = iy÷geom.ny
    ipx  = ix÷geom.nx
    rank = ipx+ipy*geom.npx+ipz*geom.npy*geom.npx
    return rank
end


#=
 Destructor for geometry data.

 @param[inout] data the geometry data structure whose storage is deallocated
=#
@inline function DeleteGeometry(geom::Geometry) 

  geom.partz_nz = nothing
  geom.partz_ids = nothing

  return
end




#endif # GEOMETRY_HPP



