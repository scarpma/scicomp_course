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

function harmonic!(x_n, dxdt_n)
  dxdt_n[1] = x_n[2]
  dxdt_n[2] = -x_n[1]
end

function harmonic_solution(t)
    x = Array{Float64,2}(undef,2,size(t,1))
    x[1,:] .= cos.(t)
    x[2,:] .= -sin.(t)
    return x
end

using BenchmarkTools

function expl_euler!(f!, x_n, x_np1, dt)
  f!(x_n, x_np1) # temporarily store dxdt_n into x_np1
  x_np1 .= x_n .+ dt .* x_np1
end

using Printf
using Plots

dim = 2
tIni = 0.
tEnd = 30.
dt = 0.02
# round and cast in julia
numTimesteps = round(Int, (tEnd - tIni) / dt) + 1

tValues = Array{Float64,1}(undef, numTimesteps)
xValues = Array{Float64,2}(undef, numTimesteps, dim)
tValues[1] = tIni
#xValues[1,:] .= [1., 1., 1.] # lorenz
xValues[1,:] .= [1., 0.]

for i=2:numTimesteps
  tValues[i] = tValues[i-1] + dt
  @views expl_euler!(harmonic!, xValues[i-1,:], xValues[i,:], dt)
  #@printf("t=%.4f: x[%d]=%.4f\n", tValues[i], i, xValues[i])
end

yValues = harmonic_solution(tValues)

l2err = sum(abs.(yValues .- xValues'))
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