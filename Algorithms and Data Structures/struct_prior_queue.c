#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
 
int *array;
int res, rep = 100, l = 0;
int extra = INT_MAX - 1;
 
struct element {
        int value;
        int index;
        int mass_index;
};
 
struct queue {
        struct element *data;
        int cap;
        int qcount;
} turn;
 
void initelem (int n, int b)
{
    int i;
    for (i = 0; i < n; i++) {
        turn.data[i].index = b;
        turn.data[i].mass_index = i;
    }
}
 
void initqueue (int size, int val)
{
    turn.data = (struct element*)malloc(size * sizeof(struct element));
    turn.cap = size;
    turn.qcount = size;
    initelem (size, val);
}
 
void swap(int a, int b)
{
        struct element elem = turn.data[a];
        turn.data[a] = turn.data[b];
        turn.data[b] = elem;
}
 
void heapify (int i, int n)
{
    int ind = 0;
    while (ind >= 0) {
        int l = 2*i + 1, r = l + 1, j = i;
        if (l < n && turn.data[l].value < turn.data[i].value)
            i = l;
        if (r < n && turn.data[r].value < turn.data[i].value)
            i = r;
        if (i == j) break;
        ind++;
        swap (i, j);
    }
}
 
int extmax (int n)
{
    //struct queue turn;
    if (0 > turn.qcount) {
        printf("devastation\n");
        return -1;
    }
    else {
        swap (0, turn.qcount);
        turn.data[turn.qcount].index++;
        heapify (0, turn.qcount);
        return 1;
    }
}
 
void inckey (int i)
{
    while (i > 0 && turn.data[(i - 1) / 2].value > turn.data[i].value) {
        swap ((i - 1) / 2, i);
        i = (i - 1) / 2;
    }
}
 
void insert(int k, int *length, int **a)
{
        array[l++] = turn.data[0].value;
        turn.qcount--;
        if (extmax (k) && turn.data[turn.qcount].index != length[turn.data[turn.qcount].mass_index]) {
                int i = turn.qcount;
                int j = a[turn.data[turn.qcount].mass_index][turn.data[turn.qcount].index];
                turn.data[turn.qcount].value = j;
                turn.qcount++;
                inckey (i);
        }
}
 
void buildheap (int n)
{
    int i = (n / 2) - 1;
    while (i >= 0) {
        heapify(i, n);
        i--;
    }
}
 
void heap_build_print (int n, int *length, int **el)
{
    //buildheap (n);
    int i = (n / 2) - 1;
    while (i >= 0) {
        heapify(i--, n);
    }
    i = n - 1;
    while (0 != turn.qcount) {
        insert (n, length, el);
        i--;
    }
    for (i = 0; i < res; i++) {
        if ((INT_MAX - 1) == array[i]) continue;
        else printf("%d ", array[i]);
    }
}
 
int main()
{
        int final_size = 0, i = 0, j, k;
        int *length, **a;
        scanf("%d", &k);
 
        initqueue (k, i);
        length = (int*)malloc(k * sizeof(int));
        a = malloc(k * sizeof(int*));
 
        for (i = 0; i < k; i++) {
                scanf("%d", &length[i]);
                if (0 == length[i]) length[i] = rep;
                final_size += length[i];
                a[i] = malloc(sizeof(int) * length[i]);
                if (i + 1 == k) res = final_size;
        }
        array = (int*)malloc(sizeof(int) * res);
       
        for (i = 0; i < k; i++) {
            if (rep != length[i]) {
                for (j = 0; j < length[i]; j++)
                scanf ("%d", &a[i][j]);
                turn.data[i].value = a[i][0];
            }
            else {
                for (j = 0; j < length[i]; j++)
                a[i][j] = extra;
                turn.data[i].value = a[i][0];
            }
        }
 
        heap_build_print(k, length, a);
 
        //mem-salvation
        for (i = 0; i < k; i++)
        free(a[i]);
        free(a);
        free(array);
        free(length);
        free(turn.data);
        return 0;
}