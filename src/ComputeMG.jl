function compute_mg!(x::Array{Float64,1}, A::HPCGSparseMatrix, r::Array{Float64,1}) #sp_coarse passed 
	A.is_mg_optimized = false
	return compute_mg_ref!(x, A , r)
end
