#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main()
{
    char i, str[1000000], symb;
    int k, j, len, alphabet['z'+1];

    scanf("%s", str);
    len = strlen(str);

    for (i = 'a'; i <= 'z'; i++) alphabet[i] = 0;

    for (k = 0; k < len; k++) {
        symb = str[k];   
        alphabet[symb]++; 
        if (symb == 0) break;
    }

    for (symb = 'a'; symb <= 'z'; symb++)
        for (j = 0; j < alphabet[symb]; j++) printf ("%c", symb);

    return 0;
}