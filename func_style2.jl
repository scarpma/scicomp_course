# function computing dxdt = f(x, t)
function f(x, dxdt, t)
  #dxdt .= [x[2], -x[1]] --> allocates!
  dxdt[1] =  x[2]
  dxdt[2] = -x[1]
  return
end

function main()
  t = 0.0
  x = zeros(2)
  dxdt = zeros(2)
  x .= [1.0, 0.0]
  @time f(x, dxdt, t)
  println(dxdt)
  return 0
end

main()