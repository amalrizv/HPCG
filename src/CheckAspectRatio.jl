#=
 @file CheckAspectRatio.cpp
 HPCG routine

using MPI
using Statistics

include("hpcg.jl")


function CheckAspectRatio(smallest_ratio, x, y, z, what, DoIo) 
  current_ratio = Statistics.min(Statistics.min(x, y), z) / (Statistics.max(Statistics.max(x, y), z))

  if current_ratio < smallest_ratio #ratio of the smallest to the largest
    if DoIo 
      HPCG_fout << "The " << what << " sizes (" << x << "," << y << "," << z <<
        ") are invalid because the ratio min(x,y,z)/max(x,y,z)=" << current_ratio <<
        " is too small (at least " << smallest_ratio << " is required)." << std::endl;
      HPCG_fout << "The shape should resemble a 3D cube. Please adjust and try again." << std::endl;
      HPCG_fout.flush();
    end

    MPI.Abort(MPI.COMM_WORLD, 127)

    return 127
  end

  return 0
end
