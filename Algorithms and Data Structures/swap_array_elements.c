#include <stdio.h>

int main()
{
    int i = 0, j = 0, sum = 0, b[j], a[i], p = 0;

    for (j = 0; j < 8; j++) {
        scanf("%d", &b[j]);
        sum = sum + b[j];
    }

    for (i = 0; i < 8; i++) {
        scanf("%d", &a[i]);
        for (j = 0; j < 8; j++) {
            if (a[i] == b[j]) {
                p++;
                b[j] = sum;
                break;
            }
        }
    }
    if (p == 8) printf("yes");
    else printf("no");
    return 0;
}