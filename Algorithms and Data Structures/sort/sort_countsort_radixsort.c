#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char src[100];
char dest[100];
int count[100] = {0};

int compare (char *str, int a, int b)
{
    int i, j;

    for (i = 0; i + a + 1 < strlen(str) && ((str[i + 1 + a] != ' ') && (str[i + 1 + a] != 9)); i++) continue;
    for (j = 0; j + b + 1 < strlen(str) && ((str[j + 1 + b] != ' ') && (str[i + 1 + a] != 9)); j++) continue;

    if (i > j) return 1;
    else return 0;
}

void restoration (int array[], int k, char *src)
{
    int i, j, m = 0;
    for (i = 0; i < k; i++) {
        for (j = 0; j + array[i] < strlen(src) && src[j + array[i]] != ' '; j++) {
            dest[m] = src[j + array[i]];
            m++;
        }
        dest[m] = ' ';
        m++;
    }
    dest[m] = 0;
    if (dest[strlen(dest) - 1] == ' ')
        dest[strlen(dest)-1] = 0;
}

int check (char *str, int i)
{
    if (dest[i] == ' ' && dest[i+1] == ' ') return 1;
    else return 0;
}

void radixsort (int array[], int k, char *src)
{
    int i,j;
    for (i = 0; i < k - 1; i++){
        for (j = i + 1; j < k; j++){
            if (compare (src, array[i], array[j]))
                count[i]++;
            else
                count[j]++;
        }
    }
}

void n_spc (char *dest)
{
    char dest2[100];
    int i, j = 0;
    for (i = 0; i < strlen(dest); i++) {
        if (dest[i] == ' ' && dest[i + 1] == ' ') continue;
        else {
            dest2[j] = dest[i + 1];
            j++;
        }
    }
    dest2[j] = 0;
    printf("%s", dest);
}

void csort (char *src, char *dest)
{
    int i, n = 0;
    int array[100];

    void word_index (char *str)
    {
        int i = 0;
        for (i = 0; i <= strlen(src); i++) {
            if ((i == 0) || ((src[i] != ' ') && (src[i - 1] == ' '))) {
                array[n] = i;
                n++;
            }
        }
    }

    word_index (src);//first chars of words

    int arr_2[n];//array of indexes of first chars

    radixsort (array, n, src);//sort

    for (i = 0; i < n; i++) arr_2[count[i]] = array[i];

    restoration(arr_2, n, src);//restoration of string
}

int main()
{
    int i = 0;
    gets (src);

    csort(src, dest);

    if (check(dest, i)) n_spc(dest);
    else         printf("%s", dest);

    return 0;
}