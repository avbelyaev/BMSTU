#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char str[22][22];

int compare(const void *a, const void *b)
{
    int i, count_a = 0, count_b = 0;
    for (i = 0; i < strlen(a); i++) if ( ((char*)a)[i] == 'a' ) count_a++;

    for (i = 0; i < strlen(b); i++) if ( ((char*)b)[i] == 'a') count_b++;

    if (count_a == count_b) return 0;
    return count_a > count_b ? 1 : -1;
}

void swap(int a, int b, void *base, int nel, size_t width)
{
    memcpy (((void*)base + width*nel + width*22), base + width*a, width);
    memcpy (base + width*a, base + width*b, width);
    memcpy (base + width*b, ((void*)base + width*nel + width*22), width);
}

void siftdown (int num_of_heap, void *base, size_t nel, size_t width)
{
    int child = num_of_heap*2, max = num_of_heap;

    if ((child < nel) && (compare(base + width*child, base + width*max) > 0))
        max = child;
    child++;
    if ((child < nel) && (compare(base + width*child, base + width*max) > 0))
        max = child;
    if (max == num_of_heap) return;

    swap(num_of_heap, max, base, nel, width);
    siftdown(max, base, nel, width);
}

void heapsort(void *base, size_t nel, size_t width, int (*compare)(const void *a, const void *b))
{
    int i = nel / 2;

    for (i; i >= 0; i--) siftdown (i, base, nel, width);
    nel--;
    for ( ; nel > 0; --nel) {
        swap (0, nel, base, nel, width);
        siftdown (0, base, nel, width);
    }
}

int main()
{
    int i, n;
    scanf("%d ", &n);
    for (i = 0; i < n; i++) gets(str[i]);

    //if (n == 1) goto exit_no_sort;

    size_t width = 22*sizeof(char);
    heapsort (str, n, width, compare);

    exit_no_sort:
    for (i = 0; i < n; i++) printf("%s\n", str[i]);

    return 0;
}