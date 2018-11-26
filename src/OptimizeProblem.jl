#= 
@file OptimizeProblem.cpp

 HPCG routine
=#

#include "OptimizeProblem.hpp"
#=
  Optimizes the data structures used for CG iteration to increase the
  performance of the benchmark version of the preconditioned CG algorithm.

  @param[inout] A      The known system matrix, also contains the MG hierarchy in attributes Ac and mgData.
  @param[inout] data   The data structure with all necessary CG vectors preallocated
  @param[inout] b      The known right hand side vector
  @param[inout] x      The solution vector to be computed in future CG iteration
  @param[inout] xexact The exact solution vector

  @return returns 0 upon success and non-zero otherwise

  @see GenerateGeometry
  @see GenerateProblem
=#
function OptimizeProblem(A, data, b, x, xexact) 

  # This function can be used to completely transform any part of the data structures.
  # Right now it does nothing, so compiling with a check for unused variables results in complaints

#if defined(HPCG_USE_MULTICOLORING)
  const nrow = A.localNumberOfRows
  colors = Vector((nrow, nrow) # value `nrow' means `uninitialized' initialized colors go from 0 to nrow-1
  totalColors = 1
  colors[1] = 0 # first point gets color 0

  # Finds colors in a greedy (a likely non-optimal) fashion.

  for i=1:nrow 
    if colors[i] == nrow # if color not assigned
      assigned = Vector(totalColors, 0)
      currentlyAssigned = 0
      const currentColIndices = A.mtxIndL[i]
      const currentNumberOfNonzeros = A.nonzerosInRow[i]

      for j=1:currentNumberOfNonzeros # scan neighbors
        curCol = currentColIndices[j]
        if curCol < i # if this point has an assigned color (points beyond `i' are unassigned)
          if assigned[colors[curCol]] == 0
            currentlyAssigned += 1
          end
          assigned[colors[curCol]] = 1 # this color has been used before by `curCol' point
        end # else  could take advantage of indices being sorted
      end # end scanning neighbours

      if currentlyAssigned < totalColors # if there is at least one color left to use
        for j=1:totalColors  # try all current colors
          if assigned[j] == 0 # if no neighbor with this color
            colors[i] = j
            break
          else 
            if colors[i] == nrow 
               colors[i] = totalColors
               totalColors += 1
	  end
        end # all colors tried
      end # no more color left to use
    end #if loop for if color not assigned 
  

    counters = Vector(totalColors)
    for i=1:nrow
      counters[colors[i]]++
    end
    old = Int64 
    old0 = Int64
    for i=2:totalColors
      old0 = counters[i]
      counters[i] = counters[i-1] + old
      old = old0
    end
    counters[1] = 0

  # translate `colors' into a permutation
    for i=1:nrow # for each color `c'
      colors[i] = counters[colors[i]]++
    end
#endif

  return 0
end

# Helper function (see OptimizeProblem.hpp for details)
function OptimizeProblemMemoryUse(A) 

  return 0.0

end
