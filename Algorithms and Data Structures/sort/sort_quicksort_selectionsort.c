#include <stdio.h>
#include <stdlib.h>

int *array;
int n;

void swap(int i, int j)
{
        int t = array[i];
        array[i] = array[j];
        array[j] = t;
}

void selectsort (int low, int high, void (*swap)(int, int))
{
    int i = low, j = high, k;
    while (j > low) {
        k = j;
        i = j - 1;
        while (i >= 0) {
            if (array[k] < array[i])
                k = i;
            i--;
        }
        swap (j, k);
        j--;
    }
}

void quicksortrec (int low, int high, int m, void (*swap)(int, int))
{//(wiki)
    if (high - low + 1 <= m) {
        selectsort(low, high, swap);
        return;
    }

    while (low != high) {

        int i = low, j = high;
        int mid = array[(low + high) / 2];

        while (i <= j) {

            while (array[i] < mid) i++;
            while (array[j] > mid) j--;

            if (i <= j) {
                swap(i, j);
                i++;
                j--;
            }

        }

        if ((j - low) < (high - i)) {
            if (low < j) quicksortrec(low, j, m, swap);
            low = i;
        }
        else {
            if (i < high) quicksortrec(i, high, m, swap);
            high = j;
        }

    }
}

void quicksort (int n, int m, void (*swap)(int, int))
{
    quicksortrec (0, n-1, m, swap);
}

int main()
{
    int i, n, m;

    scanf("%d", &n);
    scanf("%d", &m);

    array = (int*)malloc(n*sizeof(int));
    for (i = 0; i < n; i++) scanf("%d", array+i);

    quicksort (n, m, swap);
    //selectsort(0, n, swap);

    //printf("Result:\n");
    for (i = 0; i < n; i++) printf("%d ", array[i]);

    free(array);
    return 0;
}