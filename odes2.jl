using StaticArrays

#function f(x::Array{Float64,1})
#  return x.^3
#end

function lorenz!(x_n, dxdt_n)
  σ = 10
  ρ = 28
  β = 8/3
  dxdt_n[1] = σ * (x_n[2] - x_n[1])
  dxdt_n[2] = x_n[1] * (ρ - x_n[3]) - x_n[2]
  dxdt_n[3] = x_n[1] * x_n[2] - β * x_n[3]
end

function harmonic_solution(t)
    x = Array{Float64,2}(undef,2,size(t,1))
    x[1,:] .= cos.(t)
    x[2,:] .= -sin.(t)
    return x
end

function harmonic!(t, x, dxdti)
  dxdt[1] =  x[2]
  dxdt[2] = -x[1]
  return
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
    t = t + dt
    tValues[i+1] = t
    harmonic!(t, x, dxdt, i)
    @views x[:,i+1] .= x[:,i] .+ dt .* dxdt[:,i]
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
  t = t + dt
  tValues[2] = t
  harmonic!(t, x, dxdt, 1)
  @views x[:,2] .= x[:,1] .+ dt .* dxdt[:,1]

  for i=2:numTimesteps-1
    t = t + dt
    tValues[i+1] = t
    harmonic!(t, x, dxdt, i)
    @views x[:,i+1] .= x[:,i] .+ 0.5 * dt .* (3 .* dxdt[:,i] .- dxdt[:,i-1])
  end
  return (tValues, x)
end

using BenchmarkTools
using Printf
using Plots

#tValues = Array{Float64,1}(undef, numTimesteps)
#xValues = Array{Float64,2}(undef, dim, numTimesteps)
#tValues[1] = tIni
##xValues[:,1] .= [1., 1., 1.] # lorenz
#xValues[:, 1] .= [1., 0.]

#for i=2:numTimesteps
#  tValues[i] = tValues[i-1] + dt
#  expl_euler!(harmonic!, xValues[:,i-1], xValues[:,i], dt)
#  #@printf("t=%.4f: x[%d]=%.4f\n", tValues[i], i, xValues[i])
#end

tValues, xValues = @btime expl_euler(0., 1.0, [1.0, 0.0], 0.001)
yValues = harmonic_solution(tValues)
l2err = sum(abs.(yValues .- xValues))
@printf("l2 error = %.2e\n",l2err)

tValues, xValues = @btime ab2(0., 1.0, [1.0, 0.0], 0.001)
yValues = harmonic_solution(tValues)
l2err = sum(abs.(yValues .- xValues))
@printf("l2 error = %.2e\n",l2err)

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