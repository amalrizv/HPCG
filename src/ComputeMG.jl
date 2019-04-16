include("ComputeMG_ref.jl")
function compute_mg(A, r, x)
	A.is_mg_optimized = false
return compute_mg_ref(A , r , x)
end
