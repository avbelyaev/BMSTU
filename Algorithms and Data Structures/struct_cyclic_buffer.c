#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int *array, n, k;

char emp[6] = "EMPTY";
char enq[4] = "ENQ";
char deq[4] = "DEQ";

struct queue { 
    int *data;
    int head, tail;
    int qcount, cap;
} qforcopy;

struct queue InitQueue (int size)
{
    struct queue turn;
    turn.data = malloc(size * sizeof(int));
    turn.head = 0;
    turn.tail = 0;
    turn.qcount = 0;
    turn.cap = size;
    return turn;
}

void QueueEmpty(struct queue *s)
{
    if (0 == s->qcount) printf("true\n");
    else printf("false\n");
}

void EnQueue (struct queue *s, int x)
{
    struct queue turn;
    int i, prev_cap = s->cap;
    if (s->qcount == s->cap) {
        s->cap = s->cap * 2;
        s->head += prev_cap;
        s->data = realloc(s->data, prev_cap * 2 * sizeof(int));
  
        for (i = s->tail; i < prev_cap; i++) s->data[prev_cap + i] = s->data[i];

    }
    //struct queue turn;
    s->data[s->tail] = x;
    s->tail++;
    if (s->tail == s->cap) s->tail = 0;
    s->qcount++;


}

void DeQueue (struct queue *s)
{
    struct queue turn;
    int ret_val = s->data[s->head];
    s->head++;
    if (s->head == s->cap) s->head = 0;
    s->qcount--;
    printf("%d\n", ret_val);
}

int main()
{
    int enq_val, i = 4;
    char str[6];

    scanf ("%d", &n);

    struct queue turn;
    turn = InitQueue (i);

    for (i = 0; i < n; i++) {
        scanf ("%s", str);
        if (0 == strcmp(str, emp)) QueueEmpty(&turn);
        else {
            if (0 == strcmp(str, enq)) {
                scanf ("%d", &enq_val);
                EnQueue (&turn, enq_val);
            }
            else DeQueue (&turn);
        }
    }

    free(turn.data);
    free(array);
    return 0;
}