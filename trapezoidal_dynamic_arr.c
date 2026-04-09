#include <stdio.h>
#include <stdlib.h> // needed for atoi and malloc
#include <math.h>

double f(double x) {
    return pow(x, 3);
}

double trap(double a, double b, int n, double *x_values, double *f_values) {
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

int main(int argc, char *argv[]) {
    // Ensure there is at least one command line argument for the number of trapezoids
    if (argc < 2) {
        printf("Usage: %s <number_of_trapezoids>\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[1]);
    double a = 0.;
    double b = 1.;

    double *x_values = malloc(n*sizeof(double)); // array to store n points (dynamic array!)
    double *f_values = malloc(n*sizeof(double)); // array to store n points (dynamic array!)
    if (x_values == NULL || f_values == NULL) { // check if memory was allocated
        printf("Memory allocation failed\n");
        return 1;
    }

    printf("numerical integration of f(x)=x^3 from a=%.2f to b=%.2f using %d points\n",a,b,n);
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
