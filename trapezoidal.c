#include <stdio.h>
#include <math.h>

double f(double x) {
    return pow(x, 3);
}

double trap(double a, double b, int n) {
    double dx = (b-a)/n;
    double sum = 0.5 * (f(a) + f(b)); // end points contribution
    
    for (int i=1; i<n; i++) {
        double x = a + i*dx;
        sum += f(x);
    }

    return sum * dx;
}

int main() {
    int n = 1000;
    double a = 0.;
    double b = 1.;

    printf("numerical integration of f(x)=x^3 from a=%.2f to b=%.2f\n",a,b);
    double result = trap(a,b,n);
    double result_an = 0.25*(pow(b,4) - pow(a,4));
    printf("         result=%.4f\n",result);
    printf("analytic result=%.4f\n",result_an);
    printf("      abs error=%.2e\n",fabs(result - result_an));
    return 0;
}