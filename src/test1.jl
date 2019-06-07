include("test2.jl")

initial_value  =65.3
initial_arr  = zeros(5)
@show initial_value
change_from_to!(initial_value)
ires = change_arr!(initial_arr)
@show initial_value
@show ires 
#initial_value  =65.3

#change_from_to(initial_value)


