include("MixedBaseCounter.jl")

function ComputePrimeFactors(n, factors)
  d = Int64 
  sq = Int64(sqrt(n)+1)

  # remove 2 as a factor with shifts instead "/" and "%"
  while n > 1 & (n & 1) == 0  
    factors[2] = factors[2] +1
    n>>=1
  end

  # keep removing subsequent odd numbers
  for d = 3:sq 
    while 1 
      r = div(n, d)
      rr = rem(n,d)
      if rr == 0 
        factors[d] = facctors[d] +1
        n = r
        continue
      end
      break
    end
  d = d+1
  end
  if n > 1 || factors.size() == 0  # left with a prime or x==1
    factors[n] = factors[n]+1
  end
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

function ComputeOptimalShapeXYZ(xyz, x, y, z) 
  factors  = Dict{Int64, Int64} 

  ComputePrimeFactors( xyz, factors ) # factors are sorted: ascending order
  # there is at least one prime factor
  keys = collect(keys(factors))     # cache the first factor, move to the next one
  values = collect(values(factors))
  x = keys[1]
  if length(keys)>1
  y = keys[2]# try to cache the second factor in "y"
  end
  if length(factors) == 1  # only a single factor
    z = pow_i(x, factors[x] / 3)
    y = pow_i(x, factors[x] / 3 + ((factors[x] % 3) >= 2 ? 1 : 0))
    x = pow_i(x, factors[x] / 3 + ((factors[x] % 3) >= 1 ? 1 : 0))

  elseif length.factors == 2 && factors[x] == 1 && factors[y] == 1  # two distinct prime factors
    z = 1

  elseif length.factors == 2 && factors[x] + factors[y] == 3  # three prime factors, one repeated
    z = factors[x] == 2 ? x : y # test which factor is repeated

  elseif length.factors == 3 && factors[x] == 1 && factors[y] == 1 && iter->second == 1# three distinct and single prime factors
    z = keys[1]

  else  # 3 or more prime factors so try all possible 3-subsets

    distinct_factors= keys 
    count_factors =  values

    # count total number of prime factors in "c_main" and distribute some factors into "c1"
    c_main = MixedBaseCounter(count_factors, factors.size()) 
    c1 = MixedBaseCounter(count_factors, factors.size())

    # at the beginning, minimum area is the maximum area
    area= Float64 
    min_area = 2.0 * xyz + 1.0

    c1.next() 
	while ! c1.is_zero() 
     	   c2 = MixedBaseCounter(c_main, c1) # "c2" gets the factors remaining in "c_main" that "c1" doesn't have
           c2.next() 
	   while ! c2.is_zero()
        	tf1 = Int64(c1.product(distinct_factors))
        	tf2 = Int64(c2.product(distinct_factors))
        	tf3 = Int64(xyz / tf1/ tf2) # we derive the third dimension, we don't keep track of the factors it has

        	area = tf1 * Float64(tf2) + tf2 * Float(tf3) + tf1 * Float(tf3)
        	if area < min_area 
          		min_area = area
          		x = tf1
          		y = tf2
          		z = tf3
        	end
		c2.next()
      	 end
    	 c1.next()
      end
  end
end
