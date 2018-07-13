#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int *pi;

void prefix(char *str, int ind)
{
    int i, t = 0;
    pi[0] = 0;

        for (i = ind; i < strlen(str); i++) {
        //t = pi[t - 1];
	while ((t > 0) && (str[t] != str[i]))
            t = pi[t - 1];
        if (str[t] == str[i]) 
                t++;
	pi[i] = t;
	}
}

void res (char *str, char *str_2, int ind)
{
    int i, t = 0, k = strlen(str);

    for (i = ind; i < strlen(str_2); i++) {
        //t = pi[t - 1];
        while ((t > 0) && (str[t] != str_2[i]))
            t = pi[t - 1];
        if (str[t] == str_2[i]) 
                t++;
        if (t == k) {
                if (i - k + 1 >= 0) {
                      printf("%d ", i - k + 1);
                }
        }

    }
    printf("\n");
}

void kmpsubst(char *s, char *t)
{
	int i = 1, l, real_index;

	//printf("%s\n%s\n", s, t);
	char *s_res;

	l = strlen(s) + strlen(t) + 1;
	pi = (int*)malloc(l * sizeof(int));
	//s_res = (char*)malloc((l + 1) * sizeof(char));

	prefix (s, i);
        res (s, t, i-1);

        free(pi);
	//free(s_res);
}

int main(int argc, char **argv)
{
    //char a[100];
    //scanf("%s", a);
   // char b[100];
    //scanf("%s", b);
	kmpsubst(argv[1], argv[2]);
    //kmpsubst (a, b);

	return 0;
}