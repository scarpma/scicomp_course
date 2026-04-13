## lorenz system
#function f!(t, x, dxdt)
#  σ = 10
#  ρ = 28
#  β = 8/3
#  dxdt[1] = σ * (x[2] - x[1])
#  dxdt[2] = x[1] * (ρ - x[3]) - x[2]
#  dxdt[3] = x[1] * x[2] - β * x[3]
#end

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
    @views xtmp .= x[:,i] .+ 0.5*dt .* k1
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
plotly()

do_plot = true
tIni = 0.
tEnd = 30.
dt = 0.005
x0 = [1.0, 0.0]


solver = "expl_eul"
tValues, xValues = @btime expl_euler(tIni, tEnd, x0, dt)
yValues = f_solution(tValues)
l2err = sum(abs.(yValues .- xValues))
@printf("%s: l2 error = %.2e\n",solver,l2err)
if do_plot p = plot(tValues,xValues[1,:], xlabel="t", label=solver, lc=1, lw=1) end

solver = "rk2"
tValues, xValues = @btime rk2(tIni, tEnd, x0, dt)
yValues = f_solution(tValues)
l2err = sum(abs.(yValues .- xValues))
@printf("%s: l2 error = %.2e\n",solver,l2err)
if do_plot p = plot!(tValues,xValues[1,:], xlabel="t", label=solver, lc=2, lw=1) end

solver = "rk4"
tValues, xValues = @btime rk4(tIni, tEnd, x0, dt)
yValues = f_solution(tValues)
l2err = sum(abs.(yValues .- xValues))
@printf("%s: l2 error = %.2e\n",solver,l2err)
if do_plot p = plot!(tValues,xValues[1,:], xlabel="t", label=solver, lc=3, lw=1) end

solver = "ab2"
tValues, xValues = @btime ab2(tIni, tEnd, x0, dt)
yValues = f_solution(tValues)
l2err = sum(abs.(yValues .- xValues))
@printf("%s: l2 error = %.2e\n",solver,l2err)
if do_plot p = plot!(tValues,xValues[1,:], xlabel="t", label=solver, lc=4, lw=1) end



if do_plot
  p = plot!(tValues, yValues[1,:], xlabel="t", label="solution", lc=1, lw=1, ls=:dash)
  display(p)
end


## write data to file
#using DelimitedFiles
#data = hcat(tValues, xValues)  # combine into one matrix
#writedlm("ode.csv", data, ' ')

#plotly()

#p = plot( tValues, xValues[:,1], xlabel="t", ylabel="x(t)", title="Expl Euler", lc=1, lw=2)
#p = plot!(tValues, xValues[:,2], xlabel="t", ylabel="y(t)", lc=2, lw=2)
##p = plot!(tValues, xValues[:,3], xlabel="t", ylabel="z(t)", title="Explicit Euler Solution", lw=2)
#
#p = plot!(tValues, yValues[1,:], label="analytical", lc=1, ls=:dash)
#p = plot!(tValues, yValues[2,:]                    , lc=2, ls=:dash)
#display(p)

#p = plot(xValues[:,1],xValues[:,2],xValues[:,3], lw=2, xlabel="X",ylabel="Y",zlabel="Z")
#display(p)