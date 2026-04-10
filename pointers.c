#include <stdio.h>

int main() {
    int a = 10;
    int *p = &a;
    printf("p=%p, a=%d, *p=%d\n",p,a,*p);
    printf("sizeof(a)=%lu\n",sizeof(a));
    return 0;
}