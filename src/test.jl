using DataStructures

function gen_ord_set()
	o = OrderedSet{Int64}()
	a = [16, 48, 79, 560, 17, 49, 80, 561, 22, 24]
	for i = 1:10
		push!(o, a[i])
	end
	return o
end

function gen_dict(o)
	d = Dict{Int64, Int64}()
	lnr = 4096
	count = 0
	for x in o
		d[x] = lnr + count + 1
		count = count +1
	end
	return d
end

ord_set = gen_ord_set()
@show ord_set
new_dict = gen_dict(ord_set)
@show new_dict
