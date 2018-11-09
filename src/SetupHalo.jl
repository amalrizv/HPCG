 #=
@file SetupHalo.cpp

 HPCG routine
 =#
using MPI
include("SetupHalo_ref.jl")

#=
  Prepares system matrix data structure and creates data necessary necessary
  for communication of boundary values of this process.

  @param[inout] A    The known system matrix

  @see ExchangeHalo
=#
function SetupHalo(A) 

  #The call to this reference version of SetupHalo can be replaced with custom code.
  #However, any code must work for general unstructured sparse matrices.  Special knowledge about the
  #specific nature of the sparsity pattern may not be explicitly used.

  return(SetupHalo_ref(A))
end
