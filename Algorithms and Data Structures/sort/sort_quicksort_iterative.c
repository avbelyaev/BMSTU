#include <stdlib.h>
#include <stdio.h>

struct task{
    int low, high;
};

struct stack {
    struct task *array;
    int top;
    int cap;
};

void push (struct stack *st, int l, int r)
{
    st->array[st->top].low = l;
    st->array[st->top].high = r - 1;
    st->top++;
}

void pop (struct stack *stk, struct task *tsk)
{
    stk->top--;
    tsk->low = stk->array[stk->top].low;
    tsk->high = stk->array[stk->top].high;
}

void swap (int *arr, int i, int j)
{
    int temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
}

void iter_qs (int *array, int n)
{
    
    struct stack stk;
	struct task tsk;
	int i = 0, j = n - 1, mid, supp, size = n;
    stk.array = malloc(size * sizeof(struct task));
    stk.top = 0;
    stk.cap = n;
    stk.array[stk.top].low = i;
    stk.array[stk.top].high = j;
    stk.top++;

	while (stk.top != 0) {
		pop (&stk, &tsk);
		while (tsk.high > tsk.low) {
			i = tsk.low;
			j = tsk.high;
			supp = array[(i + j) / 2];
			//stdqs (array, i, j, supp);
			while (i <= j) {
                while (array[i] < supp) i++;
                while (array[j] > supp) j--;
				if (i <= j) {
					swap (array, i, j);
					i++;
					j--;
				}
			}
			if (i < ((tsk.low + tsk.high) / 2) && i < tsk.high && tsk.low >= j) {
				push (&stk, i, tsk.high + 1);
				tsk.high = j;
				continue;
			}
            push (&stk, tsk.low, j + 1);
            tsk.low = i;
		}
	}
	free (stk.array);
}

int main ()
{
    int *array;
    int i, n;

    scanf ("%d", &n);

    array = (int*)malloc(n * sizeof(int));
    for (i = 0; i < n; i++) scanf("%d", array + i);

    iter_qs (array, n);

    for (i = 0; i < n; i++) printf("%d ", array[i]);

    free(array);
    return 0;
}