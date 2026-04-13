using Printf

f(x) = x^3

#trap(a::Float64, b::Float64, n::Int) = begin
#    dx = (b - a) / n
#    x_values = a .+ (0:n-1) .* dx
#    f_values = f.(x_values)
#    s = 0.5 * (f(a) + f(b)) + sum(f_values[2:end])
#    return s * dx, x_values, f_values
#end

function trap(a::Float64, b::Float64, n::Int)
    dx = (b - a) / n
    s = 0.5 * (f(a) + f(b))

    x_values = Vector{Float64}(undef, n)
    f_values = Vector{Float64}(undef, n)

    x_values[1] = a
    f_values[1] = f(x_values[1])

    @inbounds for i in 2:n
        x_values[i] = a + (i - 1) * dx
        f_values[i] = f(x_values[i])
        s += f_values[i]
    end

    return s * dx, x_values, f_values
end

n = 1000
a, b = 0.0, 1.0

# now the function output is a tuple, not an integer
result, x_values, f_values = trap(a, b, n)
result_an = 0.25 * (b^4 - a^4)

@printf("numerical integration of f(x)=x^3 from a=%.2f to b=%.2f\n", a, b)

@inbounds for i in 1:n
    @printf("x[%4d]=%.5f, f(x[%4d])=%.5f\n", i-1, x_values[i], i-1, f_values[i])
end

@printf("         result=%.4f\n", result)
@printf("analytic result=%.4f\n", result_an)
@printf("      abs error=%.2e\n", abs(result - result_an))