include("ComputeDotProduct_ref.jl")

function compute_dot_product!(result, time_allreduce, n, x, y, is_opt)
    is_opt = false
    return compute_dot_product_ref!(result, time_allreduce, n, x, y)
end
