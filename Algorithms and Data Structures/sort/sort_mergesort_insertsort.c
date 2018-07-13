#include <stdio.h>
#include <stdlib.h>

int *array;
int *temp;
int n;

void merge(int first, int last)
{
    int i, j, k = 0, dummy, mid = (first + last)/2, first_new = first, last_new = last;

    if ((first - last + 1) < 6) {
        for (i = 1; i <= last; i++) {
            j = i;
            while ((j > 0) && ((abs(array[j])) < (abs(array[j-1])))) {
                dummy = array[j];
                array[j] = array[j - 1];
                array[j - 1] = dummy;
                j--;
            }
        }
        return;
    }

    temp = (int*)malloc(n * sizeof(int));
    merge(first, mid);
    merge(mid + 1, last);

    while (last - first + 1 != k) {
        if (first_new > mid)
            temp[k++] = array[last_new++];
        else
            if (last_new > last)
                temp[k++] = array[first_new++];
            else
                if (abs(array[first_new]) > abs(array[last_new]))
                    temp[k++] = array[last_new++];
                else
                    temp[k++] = array[first_new++];
    }

    for (i = 0; i < k; i++)
        array[first + i] = temp[i];
}

int main()
{
    int i, n;
    scanf("%d", &n);

    array = (int*)malloc(n * sizeof(int));
    for (i = 0; i < n; i++) scanf("%d", array + i);

    merge(0, n-1);

    for (i = 0; i < n; i++) printf("%d ", array[i]);

    free(array);
    free(temp);
    return 0;
}