using MPI

include("hpcg.jl")

include("ComputeSPMV.jl")
include("ComputeMG.jl")
include("ComputeDotProduct.jl")
include("ComputeResidual.jl")
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
function test_symmetry(A, b, xexact, testsymmetry_data) 

 nrow = A.localNumberOfRows
 ncol = A.localNumberOfColumns

 x_ncol = Vector{Float64}(undef, ncol)
 y_ncol = Vector{Float64}(undef, ncol)
 z_ncol = Vector{Float64}(undef, ncol)

 t4 = 0.0 # Needed for dot-product call, otherwise unused
 count_fail = 0

 # Test symmetry of matrix

 # First load vectors with random values
 fill_random_vector!(x_ncol)
 fill_random_vector!(y_ncol)

 ANorm = 2 * 26.0

 # Next, compute x'*A*y
 yNorm2, t4, ierr	= compute_dot_product!(nrow, y_ncol, y_ncol, A.is_dot_prod_optimized)
 ierr	= compute_spmv!(z_ncol, A, y_ncol) # z_nrow = A*y_overlap
 if ierr == 1 
	@debug("Error in call to SpMV: $ierr.\n")
 end

 xtAy, t4, ierr	= compute_dot_product!(nrow, x_ncol, z_ncol, A.is_dot_prod_optimized) # x'*A*y

 if ierr == 1
	@debug("Error in call to dot: $ierr .\n")
 end

 # Next, compute y'*A*x
 xNorm2, t4, ierr	= compute_dot_product!(nrow, x_ncol, x_ncol, A.is_dot_prod_optimized)
 ierr	= compute_spmv!(z_ncol, A, x_ncol) # b_computed = A*x_overlap

 if ierr == 1 
	@debug("Error in call to SpMV: $ierr .\n")
 end

 ytAx, t4, ierr	= compute_dot_product!(nrow, y_ncol, z_ncol, A.is_dot_prod_optimized) # y'*A*x

 if ierr == 1
	@debug("Error in call to dot: $err .\n")
 end
 # TODO : check if eps(Float64) is apt porting of DBL_EPSILON
 #depsym_spmv = ((xtAy - ytAx)/((xNorm2*ANorm*yNorm2 + yNorm2*ANorm*xNorm2) * eps(Float64)))

 if depsym_spmv > 1.0
	count_fail += 1  # If the difference is > 1, count it wrong
 end

 if A.geom.rank == 0
    @debug("Departure from symmetry (scaled) for SpMV abs(x'*A*y - y'*A*x) = $(depsym_spmv) ")
 end

 # Test symmetry of multi-grid

 # Compute x'*Minv*y
 ierr 	= compute_mg!(z_ncol, A, y_ncol) # z_ncol = Minv*y_ncol

 if ierr == 1
	@debug("Error in call to MG: $ierr .\n")
 end

 xtMinvy, t4, ierr		= compute_dot_product!(nrow, x_ncol, z_ncol, A.is_dot_prod_optimized) # x'*Minv*y

 if ierr == 1
	@debug("Error in call to dot: $ierr .\n")
 end

 # Next, compute z'*Minv*x
 ierr		= compute_mg!(z_ncol, A, x_ncol) # z_ncol = Minv*x_ncol

 if ierr == 1
  @debug("Error in call to MG: $ierr .\n")
 end

 ytMinvx, t4, ierr	= compute_dot_product!(nrow, y_ncol, z_ncol, A.is_dot_prod_optimized) # y'*Minv*x

 if ierr == 1 
    @debug("Error in call to dot: $ierr .\n")
 end

 depsym_mg = (xtMinvy - ytMinvx)/((xNorm2*ANorm*yNorm2 + yNorm2*ANorm*xNorm2) * (eps(Float64)))

 if depsym_mg > 1.0 
	count_fail+=1  # If the difference is > 1, count it wrong
 end

 if A.geom.rank==0
	 @debug("Departure from symmetry (scaled) for MG abs(x'*Minv*y - y'*Minv*x) = $testsymmetry_data.depsym_mg ")
 end

 x_ncol[1:length(xexact)] = xexact # Copy exact answer into overlap vector

 numberOfCalls 	= 2
 residual 	= 0.0

 for i=0:numberOfCalls 
   ierr	= compute_spmv!(z_ncol, A, x_ncol) # b_computed = A*x_overlap

   if ierr==1
	@debug("Error in call to SpMV: $ierr .\n")
   end

   residual, ierr = compute_residual!(A.localNumberOfRows, b, z_ncol)

   if ierr == 1
     @debug("Error in call to compute_residual: $ierr .\n")
   end

   if (A.geom.rank==0) 
	@debug("SpMV call [", i ,"] Residual [" ,residual, "]")
   end

 end
 testsymmetry_data = TestSymmetryData(depsym_spmv , depsym_mg , count_fail)

 x_ncol = nothing
 y_ncol = nothing
 z_ncol = nothing

 return 0
end

