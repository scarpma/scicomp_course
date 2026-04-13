#include <stdio.h>

int main() {
    int a = 10;
    int *p = &a;
    printf("p=%p, a=%d, *p=%d\n",p,a,*p);
    printf("sizeof(a)=%lu\n",sizeof(a));
    printf("as you can see, a and *p are the same thing\n");
    return 0;
}