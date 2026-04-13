# function computing dxdt = f(x, t)
function f(x, t)
  return [x[2], -x[1]]
end

function main()
  t = 0.0
  x = zeros(2)
  x .= [1.0, 0.0]
  @time dxdt = f(x, t)  # this allocates !!
  println(dxdt)
  return 0
end

main()