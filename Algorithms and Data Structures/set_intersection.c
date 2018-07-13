#include <stdio.h>
 
int main()
{
        int sizea;
        int i = 0, na, nb, resa = 0;
        scanf ("%d", &sizea);
        for (i = 0; i < sizea; i++) {
                scanf ("%d", &na);
                resa = resa | (1 << na);
        }
 
        int sizeb, resb = 0;
        scanf ("%d", &sizeb);
        for (i = 0; i < sizeb; i++) {
                scanf ("%d", &nb);
                resb = resb | (1 << nb);
        }
 
        unsigned int intersect = 0;
        intersect = resa & resb;
        i = 0;
        while (i <= 31) {
                if ((intersect & (1 << i)) > 0) printf ("%d ", i);
                i++;
        }
        return 0;
}