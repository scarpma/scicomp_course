def f(x):
    return x**3.

def trap(a,b,n):
    dx = (b-a)/n
    sum = 0.5 * (f(a)+f(b))

    for i in range(1,n):
        x = a + i*dx
        sum += f(x)

    return sum * dx

n = 1000
a = 0.
b = 1.

result = trap(a,b,n)
result_an = 0.25*(b**4-a**4)

print(f"         result={result:.4f}")
print(f"analytic result={result_an:.4f}")
print(f"      abs error={abs(result - result_an):.2e}")