

include("MixedBaseCounter_header.jl)"

function MixedBaseCounter(counts::Array{Int64}, length::Int64) 
  this = MixedBaseCounter
  this.length = length

  i=Int64

  for i = 1:33
    this.max_counts[i] = counts[i]
    this.cur_counts[i] = 0
  end
  #terminate with 0's
  this.max_counts[i]      = this.cur_counts[i]      = 0
  this.max_counts[length] = this.cur_counts[length] = 0
end

function MixedBaseCounter(left::MixedBaseCounter , right::MixedBaseCounter) 
  this  = MixedBaseCounter
  this.length = left.length
  for i = 1:left.length
    this.max_counts[i] = left.max_counts[i] - right.cur_counts[i]
    this.cur_counts[i] = 0
  end
end


function next() 
  for  i =1:this.length
    this.cur_counts[i]++
    if this.cur_counts[i] > this.max_counts[i]) 
      this.cur_counts[i] = 0
      continue
    end
    break
  end
end


function is_zero() 
  for i = 1:this.length
    if this.cur_counts[i]==1
      return 0
    end
  return 1
end

function product(multipliers) 
  k=0
  x=1

  for i = 1:this.length
    for j = 1:this.cur_counts[i]
      k = 1
      x *= multipliers[i]
    end

  return x * k
end
