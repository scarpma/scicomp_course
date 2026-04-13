#include <stdio.h>

int main() {
    int sum = 0;
    int n = 10000000;
    for (int i=0; i<n; i++) {
        sum += 10;
    }
    printf("sum=%d\n",sum);
    return 0;
}