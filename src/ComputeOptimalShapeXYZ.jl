
function compute_prime_factors(n)
  factors = Dict()
  d = Int64 
  sq = round(Int64, (sqrt(n)+1))

  # remove 2 as a factor with shifts instead "/" and "%"
    factors[2] = 0
  while n > 1 && (n & 1) == 0  
    factors[2] = factors[2] + 1
    n>>=1
  end
  # keep removing subsequent odd numbers
  for d = 3:sq 
    factors[d] = 0
    while true 
      r = div(n, d)
      rr = rem(n,d)
      if rr == 0 
        factors[d] = factors[d] +1
        n = r
        continue
      end
      break
    end
  d = d+1
  end
  factors[n] = 0
  if n > 1 || length(factors) == 0  # left with a prime or x==1
     factors[n] = factors[n]+1
  end
  for (k,v) in factors 
	if v== 0
		pop!(factors, k)
	end
  end
  return factors
end

function pow_i(x,p) 
  v = Int64

  if 0 == x || 1 == x
	 return x
  end

  if p < 0
    return 0
  end
  v = 1
  while p!=0
    if 1 & p == 1
      v *= x
    end
    x *= x
    p >>= 1
  end

  return v
end

function compute_optimal_shape_xyz(xyz, x, y, z) 
  factors = compute_prime_factors(xyz) # factors are sorted: ascending order
  if xyz == 1
	factors[1] = 1
  end
  # there is at least one prime factor
  keyss = sort!(collect(keys(factors)))     # cache the first factor, move to the next one
  vals = Array{Int64, 1}(undef, length(keyss))
  for i = 1:length(keyss)
    vals[i] = factors[keyss[i]]
  end
  x = keyss[1]

  if length(keyss)>1
      y = keyss[2] # try to cache the second factor in "y"
  end
  if length(factors) == 1  # only a single factor
    z = pow_i(x, factors[x] ÷ 3)
    y = pow_i(x, factors[x] ÷ 3 + ((factors[x] % 3) >= 2 ? 1 : 0))
    x = pow_i(x, factors[x] ÷ 3 + ((factors[x] % 3) >= 1 ? 1 : 0))
  elseif length(factors) == 2 && factors[x] == 1 && factors[y] == 1  # two distinct prime factors
    z = 1
  elseif length(factors) == 2 && factors[x] + factors[y] == 3  # three prime factors, one repeated
    z = factors[x] == 2 ? x : y # test which factor is repeated
  elseif length(factors) == 3 && factors[x] == 1 && factors[y] == 1 && keyss[3] == 1# three distinct and single prime factors
    z = keyss[1]
  else  # 3 or more prime factors so try all possible 3-subsets
    distinct_factors = Array{Int64,1}(undef, 33)
    count_factors = Array{Int64,1}(undef, 33)
    for i =1: length(keyss)
	distinct_factors[i] = keyss[i] 
    	count_factors[i] = vals[i]
    end

    # count total number of prime factors in "c_main" and distribute some factors into "c1"
    c_main = MBCounter(count_factors, length(factors))
    c1     = MBCounter(count_factors, length(factors))

    # at the beginning, minimum area is the maximum area
    area = 2.0 *xyz +1.0 
    min_area = 2.0 * xyz + 1.0

    next(c1) 

	while is_zero(c1)==0
     	   c2 = MBCounter_lr(c_main, c1) # "c2" gets the factors remaining in "c_main" that "c1" doesn't have
           next(c2) 
	   while is_zero(c2)==0

        	tf1 = Int64(product(c1,distinct_factors))
        	tf2 = Int64(product(c2,distinct_factors))
        	tf3 = Int64(xyz ÷ tf1÷ tf2) # we derive the third dimension, we don't keep track of the factors it has
        	area = tf1 * (tf2) + tf2 *(tf3) + tf1 *(tf3)
        	if area < min_area 
          		min_area = area
          		x = tf1
          		y = tf2
          		z = tf3
        	end
            	next(c2)
      	   end
    	   next(c1)
         end
  end

  return x, y, z 
end
