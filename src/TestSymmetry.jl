using MPI

include("hpcg.jl")

include("ComputeSPMV.jl")
include("ComputeMG.jl")
include("ComputeDotProduct.jl")
include("ComputeResidual.jl")
include("Geometry.jl")
include("SpMatrix.jl")
include("TestSymmetry_struct.jl")

#=
  Tests symmetry-preserving properties of the sparse matrix vector multiply and multi-grid routines.

  @param[in]    geom   The description of the problem's geometry.
  @param[in]    A      The known system matrix
  @param[in]    b      The known right hand side vector
  @param[in]    xexact The exact solution vector
  @param[inout] testsymmetry_data The data structure with the results of the CG symmetry test including pass/fail information

  @return returns 0 upon success and non-zero otherwise

  @see compute_dot_product!
  @see compute_dot_product_ref!
  @see ComputeSPMV
  @see ComputeSPMV_ref
  @see ComputeMG
  @see ComputeMG_ref
=#
function TestSymmetry(A, b, xexact, testsymmetry_data) 

 nrow = A.localNumberOfRows
 ncol = A.localNumberOfColumns

 x_ncol = Vector(undef, ncol)
 y_ncol = Vector(undef, ncol)
 z_ncol = Vector(undef, ncol)

 t4 = 0.0 # Needed for dot-product call, otherwise unused
 testsymmetry_data.count_fail = 0

 # Test symmetry of matrix

 # First load vectors with random values
 FillRandomVector(x_ncol)
 FillRandomVector(y_ncol)

 xNorm2= Float64
 yNorm2= Float64
 ANorm = 2 * 26.0

 # Next, compute x'*A*y
 compute_dot_product!(yNorm2, t4, nrow, y_ncol, y_ncol, A.isDotProductOptimized)
 ierr = ComputeSPMV(A, y_ncol, z_ncol) # z_nrow = A*y_overlap
 if ierr == 1 
	@debug("Error in call to SpMV: $ierr.\n")
 end

 xtAy = 0.0
 ierr = compute_dot_product!(xtAy, t4, nrow, x_ncol, z_ncol, A.isDotProductOptimized) # x'*A*y

 if ierr == 1
	@debug("Error in call to dot: $ierr .\n")
 end

 # Next, compute y'*A*x
 compute_dot_product!(xNorm2, t4, nrow, x_ncol, x_ncol, A.isDotProductOptimized)
 ierr = ComputeSPMV(A, x_ncol, z_ncol) # b_computed = A*x_overlap

 if ierr == 1 
	@debug("Error in call to SpMV: $ierr .\n")
 end

 ytAx = 0.0
 ierr = compute_dot_product!(ytAx, t4, nrow, y_ncol, z_ncol, A.isDotProductOptimized) # y'*A*x

 if ierr == 1
	@debug("Error in call to dot: $err .\n")
 end

 testsymmetry_data.depsym_spmv = (xtAy - ytAx)/((xNorm2*ANorm*yNorm2 + yNorm2*ANorm*xNorm2) * (DBL_EPSILON))

 if testsymmetry_data.depsym_spmv > 1.0
	testsymmetry_data.count_fail+=1  # If the difference is > 1, count it wrong
 end

 if A.geom.rank==0
    @debug("Departure from symmetry (scaled) for SpMV abs(x'*A*y - y'*A*x) = $(testsymmetry_data.depsym_spmv) ")
 end

 # Test symmetry of multi-grid

 # Compute x'*Minv*y
 ierr = ComputeMG(A, y_ncol, z_ncol) # z_ncol = Minv*y_ncol

 if ierr
	@debug("Error in call to MG: $ierr .\n")
 end

 xtMinvy = 0.0
 ierr = compute_dot_product!(xtMinvy, t4, nrow, x_ncol, z_ncol, A.isDotProductOptimized) # x'*Minv*y

 if ierr 
	@debug("Error in call to dot: $ierr .\n")
 end

 # Next, compute z'*Minv*x
 ierr = ComputeMG(A, x_ncol, z_ncol) # z_ncol = Minv*x_ncol

 if ierr 
  @debug("Error in call to MG: $ierr .\n")
 end

 ytMinvx = 0.0
 ierr = compute_dot_product!(ytMinvx, t4, nrow, y_ncol, z_ncol, A.isDotProductOptimized) # y'*Minv*x

 if ierr 
    @debug("Error in call to dot: $ierr .\n")
 end

 testsymmetry_data.depsym_mg = (xtMinvy - ytMinvx)/((xNorm2*ANorm*yNorm2 + yNorm2*ANorm*xNorm2) * (DBL_EPSILON))

 if testsymmetry_data.depsym_mg > 1.0 
	testsymmetry_data.count_fail+=1  # If the difference is > 1, count it wrong
 end

 if A.geom.rank==0
	 @debug("Departure from symmetry (scaled) for MG abs(x'*Minv*y - y'*Minv*x) = $testsymmetry_data.depsym_mg ")
 end

 xncol = xexact # Copy exact answer into overlap vector

 numberOfCalls = 2
 residual = 0.0

 for i=0:numberOfCalls 
   ierr = ComputeSPMV(A, x_ncol, z_ncol) # b_computed = A*x_overlap

   if ierr==1
	@debug("Error in call to SpMV: $ierr .\n")
   end

   ierr = ComputeResidual(A.localNumberOfRows, b, z_ncol, residual)

   if err==1
     @debug("Error in call to compute_residual: $ierr .\n")
   end

   if (A.geom.rank==0) 
	@debug("SpMV call [", i ,"] Residual [" ,residual, "]")
   end

 end

 free(x_ncol)
 free(y_ncol)
 free(z_ncol)

 return 0
end

