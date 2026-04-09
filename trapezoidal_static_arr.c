#include <stdio.h>
#include <math.h>

double f(double x) {
    return pow(x, 3);
}

double trap(double a, double b, int n, double x_values[], double f_values[]) {
    double dx = (b-a)/n;
    double sum = 0.5 * (f(a) + f(b)); // end points contribution

    x_values[0] = a;
    f_values[0] = f(x_values[0]);
    
    for (int i=1; i<n; i++) {
        double x = a + i*dx;
        x_values[i] = a + i * dx;
        f_values[i] = f(x_values[i]);  // Store the function evaluation
        sum += f_values[i];
    }

    return sum * dx;
}

int main() {
    int n = 1000;
    double a = 0.;
    double b = 1.;

    double x_values[n]; // array to store n points (static array!)
    double f_values[n]; // array to store n points (static array!)

    printf("numerical integration of f(x)=x^3 from a=%.2f to b=%.2f\n",a,b);
    double result = trap(a,b,n,x_values,f_values);
    double result_an = 0.25*(pow(b,4) - pow(a,4));
    
    for (int i = 0; i < n; i++) {
        printf("x[%4d]=%.5f, f(x[%4d])=%.5f\n", i, x_values[i], i, f_values[i]);
    }
    
    printf("         result=%.4f\n",result);
    printf("analytic result=%.4f\n",result_an);
    printf("      abs error=%.2e\n",fabs(result - result_an));
    return 0;
}
