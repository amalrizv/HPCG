#=
@file TestSymmetry.hpp

 HPCG data structures for symmetry testing
=#


include("hpcg.jl")
include("SpMatrix.jl")
include("CGData_struct.jl")

mutable struct TestSymmetryData
  depsym_spmvi::Float64;  //!< departure from symmetry for the SPMV kernel
  depsym_mg::Float64; //!< departure from symmetry for the MG kernel
   count_fail::FLoat64;   //!< number of failures in the symmetry tests
end
