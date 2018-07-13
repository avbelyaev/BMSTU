#include <stdio.h>
#define maxfib 5000
 
int main()
{
    long x;
    scanf("%ld", &x);
    if (x > 1) {
    unsigned long f[maxfib];
    int i = 2, k = 0, m;
    f[0] = 1;
    f[1] = 1;
    m = 0;
    for (i; ; i++) {
        f[i] = f[i-1] + f[i - 2];
        m++;
        if (f[i] > x)
           goto exit1;
    }
    exit1:
 
    i -= 1;
    for (i; i != 0; i--) {
        if (f[i] <= x) {
            printf("1");
            x -= f[i];
        }
        else printf("0");
    }
    }
    else printf("%ld", x);
    //exit2:
    return 0;
}