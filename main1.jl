using BenchmarkTools

f(n) = n*n + 20

n = 100000
v = Array{Float64,1}(undef, n)
@btime for i=1:n
  v[i] = f(i)
end