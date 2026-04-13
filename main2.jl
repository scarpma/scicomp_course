using BenchmarkTools

f(n) = n*n + 20

function mainLoop()
  n = 100000
  v = Array{Float64,1}(undef, n)
  for i=1:n
    v[i] = f(i)
  end
  return 0
end

@btime mainLoop()