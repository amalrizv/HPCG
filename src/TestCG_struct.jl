#= @file TestCG.hpp

 HPCG data structure
=#

include("hpcg.jl")
include("SpMatrix.jl")
include("Vector.jl")
include("CGData.jl")


mutable struct TestCGData 
  int count_pass; #!< number of succesful tests
  int count_fail;  #!< number of succesful tests
  int expected_niters_no_prec; #!< expected number of test CG iterations without preconditioning with diagonally dominant matrix (~12)
  int expected_niters_prec; #!< expected number of test CG iterations with preconditioning and with diagonally dominant matrix (~1-2)
  int niters_max_no_prec; #!< maximum number of test CG iterations without predictitioner
  int niters_max_prec; #!< maximum number of test CG iterations without predictitioner
  double normr; #!< residual norm achieved during test CG iterationse
end



