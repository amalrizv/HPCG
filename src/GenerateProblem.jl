#=
@file GenerateProblem.cpp

 HPCG routine
=#

using MPI



#=
  Routine to generate a sparse matrix, right hand side, initial guess, and exact solution.

  @param[in]  A        The generated system matrix
  @param[inout] b      The newly allocated and generated right hand side vector (if b!=0 on entry)
  @param[inout] x      The newly allocated solution vector with entries set to 0.0 (if x!=0 on entry)
  @param[inout] xexact The newly allocated solution vector with entries set to the exact solution (if the xexact!=0 non-zero on entry)

  @see GenerateGeometry
=#

function generate_problem!(A::HPCGSparseMatrix) 

  # The call to this reference version of GenerateProblem can be replaced with custom code.
  # However, the data structures must remain unchanged such that the CheckProblem function is satisfied.
  # Furthermore, any code must work for general unstructured sparse matrices.  Special knowledge about the
  # specific nature of the sparsity pattern may not be explicitly used.

  return generate_problem_ref!(A)
end
