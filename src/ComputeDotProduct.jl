
function compute_dot_product!(n::Int64, x::Array{Float64,1}, y::Array{Float64,1})
    is_opt = false
	return compute_dot_product_ref!(n, x, y)
end
