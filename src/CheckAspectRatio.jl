#=
 @file CheckAspectRatio.cpp
 HPCG routine
=#
using MPI
using Statistics

include("hpcg.jl")

function check_aspect_ratio(smallest_ratio, x, y, z, what, do_io) 

  current_ratio = Statistics.min(Statistics.min(x, y), z) / (Statistics.max(Statistics.max(x, y), z))

  if current_ratio < smallest_ratio # ratio of the smallest to the largest
    if do_io 
      @debug("The $what  sizes ($x $y $z) are invalid because the ratio min(x,y,z)/max(x,y,z)=$current_ratio is too small (at least $smallest_ratio  is required).")
      @debug("The shape should resemble a 3D cube. Please adjust and try again.")
    end

@static if MPI.Initialized()
    MPI.Abort(MPI.COMM_WORLD, 127)
end

    return 127
  end

  return 0

end
