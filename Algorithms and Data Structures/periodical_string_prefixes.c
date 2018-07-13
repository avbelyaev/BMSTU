#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
    int pi[40] = {0}, i, t = 0;
    char *str = argv[1];
    for (i = 1; i <= strlen(str); i++) {
        while ((t != 0) && (str[t] != str[i])) t = pi[t - 1];
        if (str[t] == str[i]) t++;
        pi[i] = t;
        if (i > 1 && (i % (i - pi[i - 1]) == 0) && (i != (i - pi[i - 1])))
            printf("%d %d\n", i, i/(i - pi[i - 1]));
        }
    return 0;
}