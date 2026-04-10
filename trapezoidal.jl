f(x) = x^3

function trap(a::Float64, b::Float64, n::Int)
    dx = (b-a)/n
    s = 0.5*(f(a) + f(b))

    # inbounds macro elides bounds check inside functions
    @inbounds for i in 1:n-1
        x = a + i*dx
        s += f(x)
    end
    return s * dx
end

a, b = 0.0, 1.0
n = 1000

result = trap(a, b, n)
result_an = 0.25 * (b^4 - a^4)

using Printf
@printf("         result=%.4f\n", result)
@printf("analytic result=%.4f\n", result_an)
@printf("      abs error=%.2e\n", abs(result - result_an))