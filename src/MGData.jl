#= 
@file MGData.jl

 HPCG data structure
=#


mutable struct MGData
  numberOfPresmootherSteps::Int64  # Call ComputeSYMGS this many times prior to coarsening
  numberOfPostsmootherSteps::Int64 # Call ComputeSYMGS this many times after coarsening
  f2cOperator::Array{Int64,1}      # 1D array containing the fine operator local IDs that will be injected into coarse space.
  rc::Array{Float64,1}                               # coarse grid residual vector
  xc::Array{Float64,1}                               # coarse grid solution vector
  Axf::Array{Float64,1}                              # fine grid residual vector
  init::Bool

  #=
   This is for storing optimized data structres created in OptimizeProblem and
   used inside optimized ComputeSPMV().
  =#
#  function optimizationData()# Do we need this?
end

function MGData()
	return MGData(0, 0, [], [], [], [], false)
end


#=
 Constructor for the data structure of CG vectors.

 @param[in] f2cOperator -
 @param[out] data the data structure for CG vectors that will be allocated to get it ready for use in CG iterations
 =#
@inline function init_mg_data(f2c_op::Array{Int64,1}, rc::Array{Float64,1}, xc::Array{Float64,1}, axf::Array{Float64,1}) 
  return MGData(1, 1, f2c_op, rc, xc, axf, true)
end
