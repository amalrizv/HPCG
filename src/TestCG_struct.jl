#= @file TestCG.hpp

 HPCG data structure
=#



mutable struct TestCGData 
  count_pass::Int64#!< number of succesful tests
  count_fail::Int64 #!< number of succesful tests
  expected_niters_no_prec::Int64#!< expected number of test CG iterations without preconditioning with diagonally dominant matrix (~12)
  expected_niters_prec::Int64#!< expected number of test CG iterations with preconditioning and with diagonally dominant matrix (~1-2)
  niters_max_no_prec::Int64#!< maximum number of test CG iterations without predictitioner
  niters_max_prec::Int64#!< maximum number of test CG iterations with predictitioner
  normr::Float64 #!< residual norm achieved during test CG iterationse
end



