


function MBCounter(counts::Array{Float64}, length::Int64) 
  cur_counts = zeros(33)
  this = MixedBaseCounter(length, counts,cur_counts) 
  this.length = length

  i=Int64

  for i = 1:32
    this.max_counts[i] = counts[i]
    this.cur_counts[i] = 0
  end
  i = 33
  this.cur_counts[i]   = 0
  this.max_counts[i]      = this.cur_counts[i]
  this.cur_counts[length] = 0
  this.max_counts[length] = this.cur_counts[length] 
  return this
end

function MBCounter_lr(left::MixedBaseCounter , right::MixedBaseCounter) 
  this =  MixedBaseCounter(left.length, [], [])
  temp_max = Array{Int64}(undef, left.length)
  temp_cur = Array{Int64}(undef, left.length)
  for i = 1:left.length
    temp_max[i] = left.max_counts[i] - right.cur_counts[i]
    temp_cur[i] = 0
  end
  this.max_counts = temp_max
  this.cur_counts = temp_cur
  temp_max = nothing
  temp_cur = nothing
  return this 
end


function next(this) 
  for  i =1:this.length
    this.cur_counts[i]+=1
    if this.cur_counts[i] > this.max_counts[i] 
      this.cur_counts[i] = 0
      continue
    end
    break
  end
  return this
end


function is_zero(this) 
  for i = 1:this.length
    if this.cur_counts[i]==1
      return 0
    end
   end
  return 1
end

function product(this, multipliers) 
  k=0
  x=1

  for i = 1:this.length
    for j = 1:this.cur_counts[i]
      k = 1
      x *= multipliers[i]
    end
  end
  return x * k
end
