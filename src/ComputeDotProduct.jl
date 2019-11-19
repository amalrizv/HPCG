include("ComputeDotProduct_ref.jl")

function compute_dot_product!(n, x, y, is_opt)
    is_opt = false
    return compute_dot_product_ref!(n, x, y)
end
