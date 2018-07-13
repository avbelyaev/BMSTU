#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int *pi;

void match (int ind, char *str)
{
    if (ind >= strlen(str)) {
        printf("yes");
        //return;
    }
    else {
        printf("no");
        //return;
    }
}

void prefix (char *str, int ind, int k)
{
    int i, t = 0;
    pi[0] = 0;

    for (i = ind; i < strlen(str) + 1; i++) {
        //t = pi[t - 1];
        while ((t > 0) && (str[t] != str[i]))
            t = pi[t - 1];
        if (str[t] == str[i])
                t++;
        pi[i] = t;
        if ((i >= k) && (t == 0)) break;
    }
    match(i, str);
}

void pword (char *s, char *t, int c)
{
    int i = 1;
    pi = (int*)malloc((strlen(s)+strlen(t)+1) * sizeof(int));
    //int n = c;

    prefix(s, i, c);

    /*for (i = n; i < strlen(s); i++) {
        if (pi[i] == 0) {
            printf("no");
            return;
        }
        else {
            printf("yes");
            return;
        }
    }*/


    //return;
}

int main(int argc, char **argv)
{
    int i, ls = 0;
    char *a = argv[1];
    int c = strlen(a);
    //scanf("%s", a);
    char *b = argv[2];
    //scanf("%s", b);
        //pword(argv[1], argv[2]);

        int last_symb_a = strlen(a), last_symb_b = strlen(b);
	for (i = 0; i < strlen(b); i++) {
                a[last_symb_a + i] = b[i];
                //ls++;
	}
	a[last_symb_a + last_symb_b + 1] = 0;

    pword(a, b, c);

    return 0;
}