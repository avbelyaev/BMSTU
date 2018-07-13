#include<stdio.h>
#include<stdlib.h>
#include<string.h>

char *concat(char **s, int n)
{
    int i = 0, len = 0;
    for(i = 0; i <= n; i++) len = len + strlen(s[i]);
    len = len + 1;

    char *str = (char *)malloc(len);
    strcpy(str, "");

    i = 0;
    for(i = 0; i <= n; i++) strcat(str, s[i]);
    return str;
}

int main()
{
    char *str[100];
    int i = 0, n;
    scanf ("%d", &n);

    for(i = 0; i <= n; i++) str[i] = (char *)malloc(100);

    i = 0;
    for(i = 0; i <= n; i++) gets(str[i]);

    char *temp_str = concat(str, n);
    printf("%s", temp_str);

    for(i = 0; i <= n; i++) free(str[i]);
    free(temp_str);

    return 0;
}