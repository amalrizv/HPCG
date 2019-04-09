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

function CheckProblem(AAAA::Sp_coarse, b, x, xexact) 
 # accepts sp_coarse structure. Lets break it down
  AAA = AAAA.sp_matrix		#sp_anx structure ...neighbors sendBuffer etc	
  AA = AAA.sp_matrix 		#sp_matrix structure...mIG MIL matrix values etc
  A  = AA.sp_matrix		#sp_init structure ... geom etc
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
     giz = giz0+iz
    for iy=1:ny  
       giy = giy0+iy
      for ix=1:nx  
         gix = gix0+ix
         currentLocalRow = (iz-1)*nx*ny+(iy-1)*nx+(ix-1)+1
         currentGlobalRow = giz*gnx*gny+giy*gnx+gix
        @assert(AA.localToGlobalMap[currentLocalRow] == currentGlobalRow)

        #@debug(" rank, globalRow, localRow = $A.geom.rank $currentGlobalRow ",A.globalToLocalMap[currentGlobalRow])
        numberOfNonzerosInRow = 0
        currentValuePointer = AA.matrixValues[currentLocalRow] # Pointer to current value in current row
        currentIndexPointerG = AA.mtxIndG[currentLocalRow] # Pointer to current index in current row
	cvp = 1
	cipg = 1
        for sz=-1 :1 
          if giz+sz>0 && giz+sz<=gnz
            for sy=-1:1  
              if giy+sy>0 && giy+sy<=gny 
                for sx=-1:1 
                  if gix+sx>0 && gix+sx<=gnx 
                     curcol = currentGlobalRow+sz*gnx*gny+sy*gnx+sx
                    if curcol==currentGlobalRow 
                      @assert(AA.matrixDiagonal[currentLocalRow] == currentValuePointer)
		      cvp = cvp+26.0
                     # @assert(cvp== cvp +26.0)
                     else 
		      cvp =cvp-1.0
                     # @assert(cvp== cvp -1.0)
                    end
  		      cipg  = cipg+curcol
                    #@assert(cipg== cipg+curcol)
                    numberOfNonzerosInRow+=1
                  end # end x bounds test
                end # end sx loop
              end # end y bounds test
            end # end sy loop
          end # end z bounds test
        end # end sz loop
#        @assert(AA.nonzerosInRow[currentLocalRow] == numberOfNonzerosInRow)

        localNumberOfNonzeros += numberOfNonzerosInRow # Protect this with an atomic
        if b!=0      
#		@assert(bv[currentLocalRow] == 26.0 - ((Float64)(numberOfNonzerosInRow-1)))
        end
        if x!=0      
#		@assert(xv[currentLocalRow] == 0.0)
	end
        if xexact!=0 
#		@assert(xexactv[currentLocalRow] == 1.0)
	end
      end # end ix loop
    end # end iy loop
  end # end iz loop
  @debug("Process $A.geom.rank  of $A.geom.size has $localNumberOfRows rows.\n Process $A.geom.rank of $A.geom.size has $localNumberOfNonzeros nonzeros.\n") 

   totalNumberOfNonzeros = 0
  # Use MPI's reduce function to sum all nonzeros
  #MPI.Allreduce(localNumberOfNonzeros, totalNumberOfNonzeros, MPI.SUM, MPI.COMM_WORLD)
  lnnz = localNumberOfNonzeros
  gnnz = 0 # convert to 64 bit for MPI call
  #MPI.Allreduce(lnnz, gnnz,MPI.SUM, MPI_COMM_WORLD)
  totalNumberOfNonzeros = gnnz # Copy back
  totalNumberOfNonzeros = localNumberOfNonzeros

 # @assert(AA.totalNumberOfRows == totalNumberOfRows)
 # @assert(AA.totalNumberOfNonzeros == totalNumberOfNonzeros)
 # @assert(AA.localNumberOfRows == localNumberOfRows)
 # @assert(AA.localNumberOfNonzeros == localNumberOfNonzeros)

  return
end
