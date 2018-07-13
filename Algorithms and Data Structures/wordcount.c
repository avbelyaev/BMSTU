#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int wcount (char *s)
{
    int i = 0, n = 0;
    while (i != strlen(s)){
        if (((s[i] != ' ') && (s[i] != 9))  && ((s[i+1] == ' ') || (s[i+1] == 13) || (s[i+1] == 9) || (s[i+1] == 0)))  n++;
        i++;
    }
    //printf("%d  %d  %d", s[i-1], s[i], s[i+1]);
    return n;
}

int main()
{
    char s[100];
    gets(s);
    printf ("%d", wcount (s));
   //printf("%c , %d", s[4], s[4]);
    return 0;
}