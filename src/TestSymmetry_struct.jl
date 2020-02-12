#=
@file TestSymmetry.hpp

 HPCG data structures for symmetry testing
=#



mutable struct TestSymmetryData
  depsym_spmv::Float64 #  //!< departure from symmetry for the SPMV kernel
  depsym_mg::Float64    #  //!< departure from symmetry for the MG kernel
  count_fail::Int64     #  //!< number of failures in the symmetry tests
end
