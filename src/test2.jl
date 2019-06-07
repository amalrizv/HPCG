function change_from_to!(from)
	from  = from*3-180
end
  
function change_arr!(arr)
        len = length(arr)
	for i = 1:len
		arr[i] = i*i
	end
return arr
end
