#= 
@file MGData.hpp

 HPCG data structure
=#

include("SpMatrix.jl")

mutable struct MGData
  numberOfPresmootherSteps::Int64 # Call ComputeSYMGS this many times prior to coarsening
  numberOfPostsmootherSteps::Int64 #Call ComputeSYMGS this many times after coarsening
  f2cOperator::Array{Int64,1} 1D array containing the fine operator local IDs that will be injected into coarse space.
  rc # coarse grid residual vector
  xc # coarse grid solution vector
  Axf # fine grid residual vector
  #=
   This is for storing optimized data structres created in OptimizeProblem and
   used inside optimized ComputeSPMV().
  =#
#  function optimizationData()# Do we need this?
end

#=
 Constructor for the data structure of CG vectors.

 @param[in] Ac - Fully-formed coarse matrix
 @param[in] f2cOperator -
 @param[out] data the data structure for CG vectors that will be allocated to get it ready for use in CG iterations
 =#
@inline function InitializeMGData(f2cOperator, rc, xc, Axf, data) 
  data.numberOfPresmootherSteps = 1
  data.numberOfPostsmootherSteps = 1
  data.f2cOperator = f2cOperator # Space for injection operator
  data.rc = rc
  data.xc = xc
  data.Axf = Axf
  return
end

#=
 Destructor for the CG vectors data.

 @param[inout] data the MG data structure whose storage is deallocated
=#
@inline function DeleteMGData(data) 

  free(data.f2cOperator)
  free(data.Axf)
  free(data.rc)
  free(data.xc)
  free(data.Axf)
  free(data.rc)
  free(data.xc)
  return
end


