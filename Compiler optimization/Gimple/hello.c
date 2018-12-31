#include <stdio.h>
#include <time.h>
#include <stdlib.h>
 
int notUsedFunc(int x) {
    return x * 322;
}
 
int simpleFunc(int a, int b) {
    int d = a + b;
    return d;
}
 
int main (void)
{
    srand(time(NULL));
    int cond = rand() -2;
 
    printf("Hello, World!\n");
 
    int b, c, a;
    a = 322;
    c = 633;
 
    if (cond > 0) {
        a += simpleFunc(a, b);
    } else {
        a -= simpleFunc(a, b);
    }
 
    
    printf("c=%d, a=%d\n", c, a);
    return 0;
}
