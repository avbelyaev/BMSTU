#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define max(x,y) ((x)>(y) ? (x):(y))

int fi1[1000] = {0}, fi2[1000] = {0}, sigma[500];

int *suffix (const char *str)
{
    int len = strlen(str), i = len - 2, t = i + 1;
    sigma[t] = t;

    while (i >= 0) {
        while ((t < len - 1) && (str[t] != str[i]))
            t = sigma[t + 1];
        if (str[t] == str[i]) t--;
        sigma[i] = t;
        i--;
    }
    return sigma;
}

int *delta1 (const char *str)
{
    int i = 0, j = 0;

    while (i < 26) {
        fi1[i] = strlen(str);
        i++;
    }
    while (j < strlen(str)) {
        fi1[str[j] - 61] = strlen(str) - j - 1;
        j++;
    }
    return fi1;
}

int *delta2 (const char *str)
{
    int *sigma = suffix(str);
    int len = strlen(str), i = 0, j = 0, t = sigma[0];

    while (i < len) {
        while (t < i) t = sigma[t + 1];
        fi2[i] = len + t - i;
        i++;
    }
    while (j < len - 1) {
        t = j;
        while (t < len - 1) {
            t = sigma[t + 1];
            if (str[j] != str[t]) fi2[t] = - (j + 1) + len;
        //j++;
        }
        j++;
    }
    return fi2;
}

int inc_k (const char *str_a, const char *str_b, int pos, int ind)
{
    int temp;
    if (fi1[str_b[pos] - 61] > fi2[ind]) {
        temp = pos + fi1[str_b[pos] - 61];
        return temp;
    }
    if (fi1[str_b[pos] - 61] <= fi2[ind]) {
        temp = pos + fi2[ind];
        return temp;
    }
    else return strlen(str_b);
}

int main (int argc, char **argv)
{
    //char s[100]; char t[100]; scanf("%s", s); scanf("%s", t);
    char *s = argv[1];
    char *t = argv[2];
    int tlen = strlen(t), slen = strlen(s), k = slen - 1, i;
    int *fi1 = delta1 (s);
    int *fi2 = delta2 (s);

    while (k < tlen) {
        i = slen - 1;
        cycle:;
        while (t[k] == s[i]) {
            if (i == 0) {
                printf("%d ", k);
                break;
            }
            i--;
            k--;
        }
        k = inc_k (s, t, k, i);
    }
    return 0;
}