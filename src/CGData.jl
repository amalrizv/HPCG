#=
 @file CGData.hpp
 HPCG data structure
=#


mutable struct CGData 
	r::Array{Float64,1} #< pointer to residual vector
	z::Array{Float64,1} #< pointer to preconditioned residual vector
	p::Array{Float64,1} #< pointer to direction vector
	Ap::Array{Float64,1} #< pointer to Krylov vector
end

#=
 Constructor for the data structure of CG vectors.
 @param[in]  A    the data structure that describes the problem matrix and its structure
 @param[out] data the data structure for CG vectors that will be allocated to get it ready for use in CG iterations
=#
@inline function InitializeSparseCGData(A) 
  nrow = A.localNumberOfRows
  ncol = A.localNumberOfColumns
  r    = Array{Float64,1}(undef,nrow)
  z    = Array{Float64,1}(undef,ncol)
  p    = Array{Float64,1}(undef,ncol)
  Ap   = Array{Float64,1}(undef,nrow)
  r 	= zero_fill!(r)
  z 	= zero_fill!(z)
  p 	= zero_fill!(p)
 Ap 	= zero_fill!(Ap)
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
