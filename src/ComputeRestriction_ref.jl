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
@inline function compute_restriction_ref!(A::HPCGSparseMatrix, rfv::Array{Float64,1} ) 
#  Axfv = A.mgData.Axf   #float64 array
#  rcv = A.mgData.rc     #float64 array
#  f2c = A.mgData.f2cOperator #int array
#  #KCH Performance Issues because of rcv copy
  mgd  =A.mgData
  for i::Int64  = 1:length(mgd.rc)
	@inbounds mgd.rc[i] = rfv[(mgd.f2cOperator[i])] - mgd.Axf[(mgd.f2cOperator[i])]
  end
  A.mgData = mgd
  return 0
end
