include("ComputeMG_ref.jl")
function compute_mg!(x, A, r)
	ierr = 0
	A.is_mg_optimized = false
	return compute_mg_ref!(x, A , r, ierr)
end
