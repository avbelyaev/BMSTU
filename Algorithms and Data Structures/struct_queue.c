#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define max(x,y) ((x)>(y) ? (x):(y))

int *array, n;

char emp[6] = "EMPTY";
char enq[4] = "ENQ";
char deq[4] = "DEQ";
char qmax[4] = "MAX";

struct queue { 
    int *data;
    int h;
    int q;
    int c;
    int t;
    int mc;
    int *locmax;
};

void InitQueue (struct queue *s)
{
    s->c = n;
    s->q = 0;
    s->h = n;
    s->t = 0;
    s->mc = 0;
}

void QueueEmpty(struct queue *s)
{
    if (0 != s->t) printf("false\n");
    else printf("true\n");
}

void EnQueue (struct queue *s, int value)
{
    int tmp = s->q++;
    s->data[tmp] = value;
    s->t = s->t + 1;
    if (0 == tmp || s->locmax[tmp - 1] <= value) s->locmax[tmp] = value;
    else s->locmax[tmp] = s->locmax[tmp - 1];
}

void DeQueue (struct queue *s)
{
    struct queue turn;
    int maxel;
    //q = k1; h = k2; c = n
    //s->q--;
    if (0 == s->t) {
        printf("devastation\n");
        return;
    }
    else {
        //max_search(&turn, q, n);
        if (0 != s->t && s->c == s->h) {
            int q = s->q;
            int c = s->c;
            int i = 0;

            for (i = q - 1; i >= 0; i--) {
                int j = c - q + i;
                maxel = s->data[i];
                s->data[j] = maxel;
                if (i ==
                    s->q - 1
                    || s->locmax[j + 1] <= maxel)
                    s->locmax[j] = maxel;
                else if (i != s->q - 1 || s->locmax[j + 1] >= maxel)
                    s->locmax[j] = s->locmax[j + 1];
                    else  s->q = 0;
            }
            s->q = 0;
            s->h = c - q;
        }
    }
    s->mc--;
    printf("%d\n", s->data[s->h]);
    s->h++;
    s->t--;
}

void MaxQueue (struct queue *s)
{
    int ret_val;
    if (0 != s->t && s->h == s->c)
        printf("%d\n", s->locmax[s->q - 1]);
    else {
        if (0 != s->t && 0 == s->q)
            printf("%d\n", s->locmax[s->h]);
        else {
            if (s->locmax[s->q - 1] > s->locmax[s->h])
                printf("%d\n", s->locmax[s->q - 1]);
            else printf("%d\n", s->locmax[s->h]);
        }
    }
    //printf ("%d\n", max(s->locmax[s->q - 1], s->locmax[s->h]));
}

int main()
{
    int enq_val, i;
    char str[6];

    scanf ("%d", &n);

    struct queue turn;
    turn.data = malloc(n * sizeof(int));
    turn.locmax = malloc(n * sizeof(int));

    InitQueue (&turn);

    for (i = 0; i < n; i++) {
        scanf ("%s", str);
        if (0 == strcmp(str, emp)) QueueEmpty(&turn);
        else {
            if (0 == strcmp(str, enq)) {
                scanf ("%d", &enq_val);
                EnQueue (&turn, enq_val);
            }
            else
                if (0 == strcmp(str, deq)) DeQueue (&turn);
                else MaxQueue (&turn);
        }
    }

    free(turn.data);
    free(turn.locmax);
    return 0;
}