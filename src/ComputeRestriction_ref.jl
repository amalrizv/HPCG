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
function compute_restriction_ref!(A, rfv, ierr) 
  Axfv = A.mgData.Axf   #float64 array
  rcv = A.mgData.rc     #float64 array
  f2c = A.mgData.f2cOperator #int array
  nc = length(A.mgData.rc)   #int
  #KCH Performance Issues because of rcv copy
  for i  = 1:nc
	rcv[i] = rfv[Int(f2c[i])] - Axfv[Int(f2c[i])]
  end
  A.mgData.rc = rcv
  ierr = 0
end
