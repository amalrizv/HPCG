#= @file ComputeRestriction_ref.cpp

 HPCG routine
=#

#=
  Routine to compute the coarse residual vector.

  @param[inout]  A - Sparse matrix object containing pointers to mgData->Axf, the fine grid matrix-vector product and mgData->rc the coarse residual vector.
  @param[in]    rf - Fine grid RHS.


  Note that the fine grid residual is never explicitly constructed.
  We only compute it for the fine grid points that will be injected into corresponding coarse grid points.

  @return Returns zero on success and a non-zero value otherwise.
=#
function ComputeRestriction_ref(AAAA, rf) 

  Axfv = AAAA.mgData.Axf
  rfv = rf
  rcv = AAAA.mgData.rc
  f2c = AAAA.mgData.f2cOperator
  nc = length(AAAA.mgData.rc)

  for i  = 1:nc
	#rcv[i] = rfv[f2c[i]] - Axfv[f2c[i]]
  end
  return 0
end
