#= @file GenerateProblem_ref.cpp

 HPCG routine
=#

using MPI

include("hpcg.jl")



#=
  Reference version of GenerateProblem to generate the sparse matrix, right hand side, initial guess, and exact solution.

  @param[in]  	A      The known system matrix
  		b      The newly allocated and generated right hand side vector (if b!=0 on entry)
  		x      The newly allocated solution vector with entries set to 0.0 (if x!=0 on entry)
  		xexact The newly allocated solution vector with entries set to the exact solution (if the xexact!=0 non-zero on entry)

  @see GenerateGeometry
=#

function GenerateProblem_ref(A) 

  #Make local copies of geometry information.  Use global_int_t since the RHS products in the calculations
  #below may result in global range values.
  nx = A.geom.nx
  ny = A.geom.ny
  nz = A.geom.nz
  gnx = A.geom.gnx
  gny = A.geom.gny
  gnz = A.geom.gnz
  gix0 = A.geom.gix0
  giy0 = A.geom.giy0
  giz0 = A.geom.giz0

  localNumberOfRows = nx*ny*nz # This is the size of our subblock
  # If this assert fails, it most likely means that the local_int_t is set to int and should be set to long long
  @assert(localNumberOfRows>0) # Throw an exception of the number of rows is less than zero (can happen if int overflow)
  numberOfNonzerosPerRow = 27 # We are approximating a 27-point finite element/volume/difference 3D stencil

  totalNumberOfRows = gnx*gny*gnz # Total number of grid points in mesh
  # If this assert fails, it most likely means that the global_int_t is set to int and should be set to long long
  @assert(totalNumberOfRows>0) # Throw an exception of the number of rows is less than zero (can happen if int overflow)


  # Allocate arrays that are of length localNumberOfRows
  nonzerosInRow = Array{Any}(undef,localNumberOfRows)
  mtxIndG = Array{Array{Int64,1}}(undef,localNumberOfRows)
  mtxIndL =Array{Array{Int64,1}}(undef,localNumberOfRows)
  matrixValues = Array{Array{Float64,1}}(undef,localNumberOfRows)
  matrixDiagonal = Array{Array{Float64,1}}(undef,localNumberOfRows)

  b  = Vector{Int64}(undef,localNumberOfRows)
  x =  Vector{Int64}(undef,localNumberOfRows)
  xexact = Vector{Int64}(undef,localNumberOfRows)
  bv = zeros(localNumberOfRows)
  xv = zeros(localNumberOfRows)
  xexactv = zeros(localNumberOfRows)
  if b!=0
	 bv = b # Only compute exact solution if requested
  end
  if x!=0
	 xv = x # Only compute exact solution if requested
  end
  if xexact!=0
	 xexactv = xexact # Only compute exact solution if requested
  end
  localToGlobalMap = Dict()
  globalToLocalMap = Dict()

  #Use a parallel loop to do initial assignment:
  #distributes the physical placement of arrays of pointers across the memory system
  for i=1:localNumberOfRows 
    matrixValues[i] = zeros(numberOfNonzerosPerRow)
    matrixDiagonal[i] = zeros(numberOfNonzerosPerRow)
    mtxIndG[i] = zeros(numberOfNonzerosPerRow)
    mtxIndL[i] = zeros(numberOfNonzerosPerRow)
  end
#=  if HPCG_CONTIGUOUS_ARRAYS==1	Consider making HPCG_CONTIGUOS_ARRAYS A Bool 
  for i=1:localNumberOfRows 
    mtxIndL[i] = Array{Float64}(undef,numberOfNonzerosPerRow)
  end
  for i=1:localNumberOfRows 
    matrixValues[i] = Array{Float64}(undef,numberOfNonzerosPerRow)
  end
  for i=1:localNumberOfRows
   mtxIndG[i] = Array{Float64}(undef,numberOfNonzerosPerRow)
  end
#else
=#
  for i=1:localNumberOfRows
	  mtxIndL[i] = mtxIndL[1] .+ i * numberOfNonzerosPerRow
	  matrixValues[i] = matrixValues[1] .+ i * numberOfNonzerosPerRow
	  mtxIndG[i] = mtxIndG[1] .+ i * numberOfNonzerosPerRow
  end
   localNumberOfNonzeros = 0
  # TODO:  This triply nested loop could be flattened or use nested parallelism
  for iz=1:nz
     giz = giz0+iz
     for iy=1:ny 
       giy = giy0+iy
       for  ix=1:nx 
         gix = gix0+ix
         currentLocalRow = (iz-1)*nx*ny+(iy-1)*nx+(ix-1) +1
         currentGlobalRow = giz*gnx*gny+giy*gnx+gix
         globalToLocalMap[currentGlobalRow] = currentLocalRow

         localToGlobalMap[currentLocalRow] = currentGlobalRow
         @debug(" rank, globalRow, localRow = $A.geom.rank $currentGlobalRow ",globalToLocalMap[currentGlobalRow])
         numberOfNonzerosInRow = 0
         currentValuePointer = matrixValues[currentLocalRow]  #Pointer to current value in current row
         currentIndexPointerG = mtxIndG[currentLocalRow] # Pointer to current index in current row
	 cvp = 1
	 cipg = 1
         for sz=-1:1 
          if giz+sz>0 && giz+sz<=gnz
            for sy=-1:1 
              if giy+sy>0 && giy+sy<=gny
                for sx=-1:1 
                  if gix+sx>0 && gix+sx<=gnx
                     curcol = currentGlobalRow+sz*gnx*gny+sy*gnx+sx
                    if curcol==currentGlobalRow
                      matrixDiagonal[currentLocalRow] = currentValuePointer
			
		      cvp = cvp +26
                    else 
		      cvp = cvp - 1
		      if cvp ==0
 			
		      end	
                    end
		    cipg = cipg + curcol
                    numberOfNonzerosInRow = numberOfNonzerosInRow +1
                  end #  stop x bounds test
                end # stop sx loop
              end # stop y bounds test
            end # stop sy loop
           end #stop z bounds test
        end # stop sz loop
        nonzerosInRow[currentLocalRow] = numberOfNonzerosInRow
        if b!=0
	      bv[currentLocalRow] = 26.0 - (numberOfNonzerosInRow-1)
	end
        if x!=0      
		xv[currentLocalRow] = 0.0
	end
        if xexact!=0 
		xexactv[currentLocalRow] = 1.0
	end
      end #  stop ix loop
    end # stop iy loop
  end # stop iz loop
  @debug("Process $A.geom.rank of $A.geom.size has $localNumberOfRows rows.\n Process $A.geom.rank of $A.geom.size has $localNumberOfNonzeros nonzeros.\n") 

  totalNumberOfNonzeros = 0
  # Use MPI's reduce function to sum all nonzeros
#  if USE_MPI == true
#	  MPI.Allreduce(localNumberOfNonzeros, totalNumberOfNonzeros, MPI.SUM, MPI.COMM_WORLD)
#  end
  lnnz = localNumberOfNonzeros 
  gnnz = 0 # convert to 64 bit for MPI call
#  if USE_MPI==true
#  	MPI.Allreduce(lnnz, gnnz, MPI.SUM, MPI.COMM_WORLD)
#  end
  totalNumberOfNonzeros = gnnz # Copy back
  totalNumberOfNonzeros = localNumberOfNonzeros
  # If this assert fails, it most likely means that the global_int_t is set to int and should be set to long long
  # This assert is usually the first to fail as problem size increases beyond the 32-bit integer range.
#  @assert(totalNumberOfNonzeros>0) # Throw an exception of the number of nonzeros is less than zero (can happen if int overflow)
  AA= SpMatrix(A, "0", totalNumberOfRows, totalNumberOfNonzeros, localNumberOfRows,localNumberOfRows, localNumberOfNonzeros, nonzerosInRow, mtxIndG, mtxIndL, matrixValues, matrixDiagonal, localToGlobalMap, globalToLocalMap) 
  return AA
end
