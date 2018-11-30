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
function ComputeProlongation_ref(Af,xf) 

  xfv = xf
  xcv = Af.mgData.xc
  f2c = Af.mgData.f2cOperator
  nc = Af.mgData.rc.localLength

# TODO: Somehow note that this loop can be safely vectorized since f2c has no repeated indices
  for i=1:nc
	xfv[f2c[i]] += xcv[i] # This loop is safe to vectorize
  end
  return 0
end
