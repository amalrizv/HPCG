#=
@file ComputeProlongation_ref.cpp

 HPCG routine
=#


#=
  Routine to compute the coarse residual vector.

  @param[in]  Af - Fine grid sparse matrix object containing pointers to current coarse grid correction and the f2c operator.
  @param[inout] xf - Fine grid solution vector, update with coarse grid correction.

  Note that the fine grid residual is never explicitly constructed.
  We only compute it for the fine grid points that will be injected into corresponding coarse grid points.

  @return Returns zero on success and a non-zero value otherwise.
=#
function compute_prolongation_ref!(xf::Array{Float64,1}, Af::HPCGSparseMatrix) 

	mgd = Af.mgData
#  xcv = Af.mgData.xc
#  f2c = Af.mgData.f2cOperator
#  nc = length(Af.mgData.rc)

# TODO: Somehow note that this loop can be safely vectorized since f2c has no repeated indices
for i=1:length(mgd.rc)
	@fastmath @inbounds	xf[mgd.f2cOperator[i]] += mgd.xc[i] # This loop is safe to vectorize
  end
  
  return 0 
end
