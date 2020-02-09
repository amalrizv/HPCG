
mutable struct alt_set
	size::Int64
	alt::Dict{}
end

function alt_set(sz)
	return alt_set(sz, Dict())
end

function alt_set_add!(a::alt_set, to_add::Int64)

	if haskey(a.alt,to_add) == false
		val  = a.size +1
		a.alt[to_add] = val
		a.size  = val
	end
		
end

# push places the set value at the top of the set stack
# which is why they are all placed topsy turvy.
# We know from the C version that the first entry in th
# receiveList mapping should be 16, 48, 528
