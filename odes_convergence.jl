# harmonic oscillator
function f!(t, x, dxdt)
  dxdt[1] =  x[2]
  dxdt[2] = -x[1]
  return
end
function f_solution(t)
    x = Array{Float64,2}(undef,2,size(t,1))
    x[1,:] .= cos.(t)
    x[2,:] .= -sin.(t)
    return x
end


function expl_euler(tIni, tEnd, x0, dt)
  numTimesteps = round(Int, (tEnd - tIni) / dt)
  tValues = Array{Float64,1}(undef,numTimesteps)
  x       = Array{Float64,2}(undef,size(x0,1),numTimesteps)
  dxdt    = Array{Float64,2}(undef,size(x0,1),numTimesteps)
  t = tIni
  tValues[1] = t
  x[:,1] .= x0
  for i=1:numTimesteps-1
    @views f!(t, x[:,i], dxdt[:,i])
    @views x[:,i+1] .= x[:,i] .+ dt .* dxdt[:,i]
    t = t + dt
    tValues[i+1] = t
  end
  return (tValues, x)
end

function rk2(tIni, tEnd, x0, dt)
  numTimesteps = round(Int, (tEnd - tIni) / dt)
  tValues = Array{Float64,1}(undef,numTimesteps)
  x       = Array{Float64,2}(undef,size(x0,1),numTimesteps)
  dxdt    = Array{Float64,2}(undef,size(x0,1),numTimesteps)
  k1      = Array{Float64,1}(undef,size(x0,1))
  k2      = Array{Float64,1}(undef,size(x0,1))
  xtmp    = Array{Float64,1}(undef,size(x0,1))
  t = tIni
  tValues[1] = t
  x[:,1] .= x0
  for i=1:numTimesteps-1
    # k1
    @views f!(t, x[:,i], k1)
    @views xtmp .= x[:,i] .+ dt .* k1
    # k2
    @views f!(t+dt, xtmp, k2)
    @views x[:,i+1] .= x[:,i] .+ 0.5*dt .* (k1 .+ k2)
    t = t + dt
    tValues[i+1] = t
  end
  return (tValues, x)
end

function rk4(tIni, tEnd, x0, dt)
  numTimesteps = round(Int, (tEnd - tIni) / dt)
  tValues = Array{Float64,1}(undef,numTimesteps)
  x       = Array{Float64,2}(undef,size(x0,1),numTimesteps)
  dxdt    = Array{Float64,2}(undef,size(x0,1),numTimesteps)
  k1      = Array{Float64,1}(undef,size(x0,1))
  k2      = Array{Float64,1}(undef,size(x0,1))
  k3      = Array{Float64,1}(undef,size(x0,1))
  k4      = Array{Float64,1}(undef,size(x0,1))
  xtmp    = Array{Float64,1}(undef,size(x0,1))
  t = tIni
  tValues[1] = t
  x[:,1] .= x0
  for i=1:numTimesteps-1
    # k1
    @views f!(t, x[:,i], k1)
    @views xtmp .= x[:,i] .+ 0.5 .* dt .* k1
    # k2
    f!(t + 0.5*dt, xtmp, k2)
    @views xtmp .= x[:,i] .+ 0.5 .* dt .* k2
    # k3
    f!(t + 0.5*dt, xtmp, k3)
    @views xtmp .= x[:,i] .+ dt .* k3
    # k4
    f!(t + dt, xtmp, k4)
    @views x[:,i+1] .= x[:,i] .+ (dt/6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
    t = t + dt
    tValues[i+1] = t
  end
  return (tValues, x)
end

function ab2(tIni, tEnd, x0, dt)
  numTimesteps = round(Int, (tEnd - tIni) / dt)
  tValues = Array{Float64,1}(undef,numTimesteps)
  x       = Array{Float64,2}(undef,size(x0,1),numTimesteps)
  dxdt    = Array{Float64,2}(undef,size(x0,1),numTimesteps)
  t = tIni
  tValues[1] = t
  x[:,1] .= x0

  # start the integration using explicit euler
  @views f!(t, x[:,1], dxdt[:,1])
  @views x[:,2] .= x[:,1] .+ dt .* dxdt[:,1]
  t = t + dt
  tValues[2] = t

  for i=2:numTimesteps-1
    @views f!(t, x[:,i], dxdt[:,i])
    @views x[:,i+1] .= x[:,i] .+ 0.5 * dt .* (3 .* dxdt[:,i] .- dxdt[:,i-1])
    t = t + dt
    tValues[i+1] = t
  end
  return (tValues, x)
end

using BenchmarkTools
using Printf
using Plots
#plotly()

do_plot = true
tIni = 0.
tEnd = 3.
dt = 0.005
x0 = [1.0, 0.0]

dts = 0:8
dts = 0.5 .* 0.25 .^ dts
println(dts)

p = plot()
solvers = [("expl_eul", expl_euler, 1), 
           ("ab2", ab2, 2),
           ("rk2", rk2, 2),
           ("rk4", rk4, 4)]
for (i, (solver_name, solver, order)) in enumerate(solvers)
  errs = Float64[]
  for dt in dts
    tValues, xValues = solver(tIni, tEnd, x0, dt)
    yValues = f_solution(tValues)
    maxerr = maximum(abs.(yValues .- xValues))
    push!(errs, maxerr)
  end
  ref = errs[1] .* (dts ./ dts[1]) .^ order
  plot!(dts, ref, ls=:dash, lc=i, label="", lw=3)
  plot!(dts, errs, marker=:circle, label=solver_name, lc=:black, mc=i)

end
plot!(xaxis=(:log10,(1e-6,1e-0)), xflip=true)
plot!(yaxis=:log10, xlabel="dt")
plot!(yaxis=:log10, ylabel="max error")
plot!(legend=:outertopright)
display(p)