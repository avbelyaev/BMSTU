#define _USE_MATH_DEFINES
#include <math.h>
#include <locale.h>
#include <stdio.h>

FILE *file;

double f( double x ) 
{
    return 
}

double F( double x ) 
{
    return 
}

double err = 0.0;

double Integral( double Left, double Right, long N, double (*func)(double) )
{
    int i;
    double x, dx, res = 0.0;
    dx = (Right - Left) / N;
    x = Left + dx;

    for ( i = 0; i < N; i++) {
        x = Left + (i + 1) * dx;
        double y = (func(x) + func(x + dx)) / 2 - err;
        double t = res + y;
        err = (t - res) - y;
        res = t;
    }
    res *= dx;
    return res;
}


int main ()
{
    long n;
    double L = 0.0, R = 0.0;
    file = fopen ("archres_KAHAN.csv", "w");

    double V, V0 = F(R) - F(L);
    setlocale ( LC_ALL, "");
    fprintf (file, "Num of steps;Relative Mistake;Evaluation of Mistake\n");

    for (n = 1; n < 100; n += n/50 + 1) {
        V = Integral( L, R, n, f );
        fprintf (file, "%ld;=%.15G;=%.15G\n", n, (V-V0)/V0,
               // (((R-L)*(R-L)*(R-L)*(pow(10, 10)))/(12.*n*n))/V0 );

    }

    fclose(file);
    return 0;
}
//for ( x = Left, i = 0; i < N; x+= dx, i++) {
//for ( x = Left; x < Right; x += dx) {