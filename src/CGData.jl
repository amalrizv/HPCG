#=
 @file CGData.hpp
 HPCG data structure
=#


include("SpMatrix.jl")

mutable struct CGData 
  r::Vector #< pointer to residual vector
  z::Vector #< pointer to preconditioned residual vector
  p::Vector #< pointer to direction vector
  Ap::Vector #< pointer to Krylov vector
end

#=
 Constructor for the data structure of CG vectors.
 @param[in]  A    the data structure that describes the problem matrix and its structure
 @param[out] data the data structure for CG vectors that will be allocated to get it ready for use in CG iterations
=#
@inline function InitializeSparseCGData(A, data) 
  nrow = A.localNumberOfRows
  ncol = A.localNumberOfColumns
  r = Vector{Float64}(undef,nrow)
  z = Vector{Float64}(undef,ncol)
  p = Vector{Float64}(undef,ncol)
  Ap = Vector{Float64}(undef,nrow)
  data = CGData(r,z,p, Ap)
  return data
end

#=
 Destructor for the CG vectors data.
 @param[inout] data the CG vectors data structure whose storage is deallocated
=#
@inline function DeleteCGData(data) 

  free(data.r)
  free(data.z)
  free(data.p)
  free(data.Ap)
  return
end

#endif #CGDATA_HPP
