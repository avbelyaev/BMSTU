#define _USE_MATH_DEFINES

#include <math.h>
#include <locale.h>
#include <stdio.h>

double f( double x ) {     
    return //Function to integrate
}

double F( double x ) {     
    return //AntiDerivative
}



double Integral( double Left, double Right, long N, double (*func)(double) )
{
    double x, dx, res = 0.0;
    dx = (Right - Left) / N;

    for ( x = Left; x < Right; x += dx) {
        
	res += (func(x) + func(x + dx)) / 2;     //Trapezoidal rule
    }
    res *= dx;
    return res;
}



int main ()
{
    long n;
    double L = , R = ;            //Area
    double V, V0 = F(R) – F(L);       //Precise result (Newton - Leibniz)
    setlocale ( LC_ALL, "");
    fprintf (file, "Num of steps;Relative Mistake;Evaluation of Mistake\n");

    for (n = 1; n < 100; n += n/50 + 1) {    //Num of steps
        V = Integral( L, R, n, f );               
        
	fprintf (file, "%ld;=%.15G;=%.15G\n", n, (V-V0)/V0,               //Relative Error (15 digits for double)
        //        (((R-L)*(R-L)*(R-L)*(pow(10, 10)))/(12.*n*n))/V0 );     //Error of trapezoidal rule
    }
    return 0;
}