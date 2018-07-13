#include <stdio.h>
int main()
{
    int n = 0, i = 0;
    scanf("%d", &n);
    int a[n];
    for (i = 0; i < n; i++) 
        scanf("%d", &a[i]);
 
    int k = 0;
    scanf("%d", &k);       
 
    int maxel=0;
    maxel = 0;
    for (i = 0; i < k; i++)
        maxel = maxel + a[i]; 
    //printf("%d\n", maxel);
    int maxt = maxel;

 
    for (i = k; i <= n-1; i++) {
       // printf("d= %d, i= %d,", d, i);
        maxt = maxt + a[i];  
        maxt = maxt - a[i - k]; 
       //printf(" maxt =%d\n", maxt);
        if (maxt > maxel)
            maxel = maxt;
    }
    printf("%d", maxel);  
    return 0;
}