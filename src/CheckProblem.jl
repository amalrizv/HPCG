#=
 @file CheckProblem.cpp

 HPCG routine
=#

using MPI
include("hpcg.jl")


#=
  Check the contents of the generated sparse matrix to see if values match expected contents.

  @param[in]  A      The known system matrix
  @param[inout] b      The newly allocated and generated right hand side vector (if b!=0 on entry)
  @param[inout] x      The newly allocated solution vector with entries set to 0.0 (if x!=0 on entry)
  @param[inout] xexact The newly allocated solution vector with entries set to the exact solution (if the xexact!=0 non-zero on entry)

  @see GenerateGeometry
=#

function check_problem(A::HPCGSparseMatrix, b, x, xexact) 
  # Make local copies of geometry information.  Use global_int_t since the RHS products in the calculations
  # below may result in global range values.
  nx = A.geom.nx
  ny = A.geom.ny
  nz = A.geom.nz
  gnx = A.geom.gnx
  gny = A.geom.gny
  gnz = A.geom.gnz
  gix0 = A.geom.gix0
  giy0 = A.geom.giy0
  giz0 = A.geom.giz0

  localNumberOfRows = nx*ny*nz  #This is the size of our subblock
  totalNumberOfRows = gnx*gny*gnz # Total number of grid points in mesh

  bv = Vector{Float64}
  xv = Vector{Float64}
  xexactv = Vector{Float64}
  if b!=0 
	bv = b # Only compute exact solution if requested
  end
  if x!=0 
	xv = x # Only compute exact solution if requested
  end
  if xexact!=0
	 xexactv = xexact # Only compute exact solution if requested
  end
   localNumberOfNonzeros = 0
  # TODO:  This triply nested loop could be flattened or use nested parallelism
  for iz=1:nz 
     giz = giz0+(iz-1)
    for iy=1:ny  
       giy = giy0+(iy-1)
      for ix=1:nx  
         gix = gix0+(ix-1)
         currentLocalRow = (iz-1)*nx*ny+(iy-1)*nx+(ix-1)+1
         currentGlobalRow = giz*gnx*gny+giy*gnx+gix+1
        @assert(A.localToGlobalMap[currentLocalRow] == currentGlobalRow)

		@debug(" rank, globalRow, localRow = $(A.geom.rank) $currentGlobalRow ,$(A.globalToLocalMap[currentGlobalRow])")
        numberOfNonzerosInRow = 0
        currentValuePointer = 1 # Pointer to current value in current row
        currentIndexPointerG = 1 # Pointer to current index in current row
        for sz=-1 :1 
          if giz+sz>-1 && giz+sz<gnz
            for sy=-1:1  
              if giy+sy>-1 && giy+sy<gny 
                for sx=-1:1 
                  if gix+sx>-1 && gix+sx<gnx 
                     curcol = (currentGlobalRow-1)+sz*gnx*gny+sy*gnx+sx+1
                    if curcol==currentGlobalRow 
                      @assert(A.matrixValues[currentLocalRow, currentValuePointer] == 26.0)
                     else 
		      @assert(A.matrixValues[currentLocalRow, currentValuePointer] == -1.0)
                    end
		    @assert(A.mtxIndG[currentLocalRow, currentIndexPointerG] == curcol)
                    currentValuePointer +=1
		    currentIndexPointerG +=1
                    numberOfNonzerosInRow+=1
                  end # end x bounds test
                end # end sx loop
              end # end y bounds test
            end # end sy loop
          end # end z bounds test
        end # end sz loop
        @assert(A.nonzerosInRow[currentLocalRow] == numberOfNonzerosInRow)

        localNumberOfNonzeros += numberOfNonzerosInRow # Protect this with an atomic
        if b!=0      
		@assert(bv[currentLocalRow] == 26.0 - ((Float64)(numberOfNonzerosInRow-1)))
        end
        if x!=0      
		@assert(xv[currentLocalRow] == 0.0)
	end
        if xexact!=0 
		@assert(xexactv[currentLocalRow] == 1.0)
	end
      end # end ix loop
    end # end iy loop
  end # end iz loop
  @debug("Process $(A.geom.rank)  of $(A.geom.size) has $localNumberOfRows rows.\n Process $(A.geom.rank) of $(A.geom.size) has $localNumberOfNonzeros nonzeros.\n") 

  totalNumberOfNonzeros = 0
  if MPI.Initialized()== true 
     # Use MPI's reduce function to sum all nonzeros

     totalNumberOfNonzeros = MPI.Allreduce(localNumberOfNonzeros, MPI.SUM, MPI.COMM_WORLD)
  else 
     totalNumberOfNonzeros = localNumberOfNonzeros
  end 
 # @assert(A.totalNumberOfRows == totalNumberOfRows)
 # @assert(A.totalNumberOfNonzeros == totalNumberOfNonzeros)
 # @assert(A.localNumberOfRows == localNumberOfRows)
 # @assert(A.localNumberOfNonzeros == localNumberOfNonzeros)

  return
end
