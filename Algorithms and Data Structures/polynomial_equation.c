#include <stdio.h>

int main()
{
    int i = 0, n;
    long int x, a, p, pr;
    scanf ("%d%ld%ld", &n, &x, &a);
    p = a; pr = a*n;
    for (i = n; i > 0; i--) {
        scanf ("%ld", &a);
        p = p*x + a;
        if (i > 1) pr = a*(i - 1) + pr*x;
    }
    printf ("%ld %ld", p, pr);
    return 0;
}