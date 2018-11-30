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

function CheckProblem(A, b, x, xexact) 

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

  bv = 0
  xv = 0
  xexactv = 0
  if b!=0 
	bv = b.values # Only compute exact solution if requested
  end
  if x!=0 
	xv = x.values # Only compute exact solution if requested
  end
  if xexact!=0
	 xexactv = xexact.values # Only compute exact solution if requested
  end
   localNumberOfNonzeros = 0
  # TODO:  This triply nested loop could be flattened or use nested parallelism
  for iz=0 iz<nz iz++) 
     giz = giz0+iz
    for ( iy=0 iy<ny iy++) 
       giy = giy0+iy
      for ( ix=0 ix<nx ix++) 
         gix = gix0+ix
         currentLocalRow = iz*nx*ny+iy*nx+ix
         currentGlobalRow = giz*gnx*gny+giy*gnx+gix
        @assert(A.localToGlobalMap[currentLocalRow] == currentGlobalRow)

        @debug(" rank, globalRow, localRow = $A.geom.rank $currentGlobalRow  $(A.globalToLocalMap[$currentGlobalRow])")  

        numberOfNonzerosInRow = 0
        currentValuePointer = A.matrixValues[currentLocalRow] # Pointer to current value in current row
        currentIndexPointerG = A.mtxIndG[currentLocalRow] # Pointer to current index in current row
        for (int sz=-1 sz<=1 sz++) 
          if (giz+sz>-1 && giz+sz<gnz) 
            for (int sy=-1 sy<=1 sy++) 
              if (giy+sy>-1 && giy+sy<gny) 
                for (int sx=-1 sx<=1 sx++) 
                  if (gix+sx>-1 && gix+sx<gnx) 
                     curcol = currentGlobalRow+sz*gnx*gny+sy*gnx+sx
                    if (curcol==currentGlobalRow) 
                      @assert(A.matrixDiagonal[currentLocalRow] == currentValuePointer)
                      @assert(*currentValuePointer++ == 26.0)
                     else 
                      @assert(*currentValuePointer++ == -1.0)
                    end
                    @assert(*currentIndexPointerG++ == curcol)
                    numberOfNonzerosInRow++
                  end # end x bounds test
                end # end sx loop
              end # end y bounds test
            end # end sy loop
          end # end z bounds test
        end # end sz loop
        @assert(A.nonzerosInRow[currentLocalRow] == numberOfNonzerosInRow)

        localNumberOfNonzeros += numberOfNonzerosInRow # Protect this with an atomic
        if (b!=0)      @assert(bv[currentLocalRow] == 26.0 - ((double) (numberOfNonzerosInRow-1)))
        if (x!=0)      @assert(xv[currentLocalRow] == 0.0)
        if (xexact!=0) @assert(xexactv[currentLocalRow] == 1.0)
      end # end ix loop
    end # end iy loop
  end # end iz loop
  @debug("Process $A.geom.rank  of $A.geom.size has $localNumberOfRows rows.\n Process $A.geom.rank of $A.geom.size has $localNumberOfNonzeros nonzeros.\n" 

   totalNumberOfNonzeros = 0
  # Use MPI's reduce function to sum all nonzeros
  MPI.Allreduce(localNumberOfNonzeros, totalNumberOfNonzeros, MPI.SUM, MPI.COMM_WORLD)
  lnnz = localNumberOfNonzeros, gnnz = 0 # convert to 64 bit for MPI call
  MPI.Allreduce(lnnz, gnnz,MPI.SUM, MPI_COMM_WORLD)
  totalNumberOfNonzeros = gnnz # Copy back
  totalNumberOfNonzeros = localNumberOfNonzeros

  @assert(A.totalNumberOfRows == totalNumberOfRows)
  @assert(A.totalNumberOfNonzeros == totalNumberOfNonzeros)
  @assert(A.localNumberOfRows == localNumberOfRows)
  @assert(A.localNumberOfNonzeros == localNumberOfNonzeros)

  return
end
