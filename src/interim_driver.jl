using Distributed
using MPI
using Revise 
ENV["JULIA_NUM_THREADS"] = 1
Distributed.@everywhere includet("header.jl")

#Struct Procedures

Distributed.@everywhere includet("appendx.jl")
Distributed.@everywhere includet("Geometry.jl")
Distributed.@everywhere includet("hpcg.jl")
Distributed.@everywhere includet("CGData.jl")
Distributed.@everywhere includet("MGData.jl")
Distributed.@everywhere includet("SpMatrix.jl")
Distributed.@everywhere includet("TestCG_struct.jl")
Distributed.@everywhere includet("TestSymmetry_struct.jl")
Distributed.@everywhere includet("TestNorms_struct.jl")
Distributed.@everywhere includet("MixedBaseCounter_header.jl")

#Leaf procedures

Distributed.@everywhere includet("init.jl")
Distributed.@everywhere includet("MixedBaseCounter.jl")
Distributed.@everywhere includet("ExchangeHalo.jl")
Distributed.@everywhere includet("ComputeDotProduct_ref.jl")
#Distributed.@everywhere includet("ComputeDotProduct.jl")
Distributed.@everywhere includet("ComputeOptimalShapeXYZ.jl")
Distributed.@everywhere includet("ComputeProlongation_ref.jl")
Distributed.@everywhere includet("ComputeResidual.jl")
Distributed.@everywhere includet("ComputeRestriction_ref.jl")
Distributed.@everywhere includet("ComputeRestriction.jl")
Distributed.@everywhere includet("ComputeSYMGS_ref.jl")
Distributed.@everywhere includet("ComputeWAXPBY_ref.jl")
Distributed.@everywhere includet("ComputeWAXPBY.jl")

#Procedures in  main

Distributed.@everywhere includet("GenerateGeometry.jl")
Distributed.@everywhere includet("CheckAspectRatio.jl")
Distributed.@everywhere includet("GenerateProblem_ref.jl")
Distributed.@everywhere includet("GenerateProblem.jl")
Distributed.@everywhere includet("SetupHalo_ref.jl")
Distributed.@everywhere includet("SetupHalo.jl")
Distributed.@everywhere includet("GenerateCoarseProblem.jl")
Distributed.@everywhere includet("CheckProblem.jl")
Distributed.@everywhere includet("ComputeSPMV_ref.jl")
Distributed.@everywhere includet("ComputeSPMV.jl")
Distributed.@everywhere includet("ComputeMG_ref.jl")
Distributed.@everywhere includet("ComputeMG.jl")
Distributed.@everywhere includet("OptimizeProblem.jl")
#includet("mytimer.jl")
Distributed.@everywhere includet("CG_ref.jl")
Distributed.@everywhere includet("CG.jl")
Distributed.@everywhere includet("TestCG.jl")
Distributed.@everywhere includet("TestSymmetry.jl")
Distributed.@everywhere includet("TestNorms.jl")
Distributed.@everywhere includet("ReportResults.jl")

#main

Distributed.@everywhere includet("main.jl")

hpcg_args = header_calling_hpcg()
precompile(compute_dot_product_ref!, (Int64, Array{Float64,1}, Array{Float64,1}, Float64, Float64))
precompile(compute_spmv_ref!, (Array{Float64,1}, HPCGSparseMatrix, Array{Float64,1}))
precompile(compute_mg_ref!, (Array{Float64,1}, HPCGSparseMatrix, Array{Float64,1})) 
# KCH: DEBUGGING SHOULD NOT BE ON FOR PERFORMANCE TESTING
#ENV["JULIA_DEBUG"] = "all" 
main(hpcg_args)

