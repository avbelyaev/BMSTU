#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
 
int n, comp_sign, extra = INT_MAX;

struct Elem {
  struct Elem *prev, *next;
  int v;
};

struct Elem *x, *t1, *el, *z, *dy, *element, *current, *dz, *l, *check, *start, *finish, *res, *elem, *temp1, *temp2, *print;

struct Elem* initlist()
{
  struct Elem *first=(struct Elem*)malloc(sizeof(struct Elem));
  first->prev=first;
  first->next=first;
  return first;
}

struct Elem* listsearch (struct Elem *ls, int v)
{
    x = ls->next;
    while (x != ls && x->v != v)
        x = x->next;
    return x;
}

struct Elem* listlength (struct Elem *ls) {
    int len = 0;
    x = ls;
    while (x->next != ls) {
        len++;
        x = x->next;
    }
    if (len != n) printf("loss\n");
    return len;
}

void Insert(struct Elem *a, struct Elem *b){
  b->next=a->next;
  b->prev=a;
  a->next=b;
  b->next->prev=b;
}

int compare (struct Elem *x, struct Elem *y)
{
    if (x->v < y->v) {
        //comp_sign = 1;
        return 1;
    }
    else return 0;
    //comp_sign = 0;
    //return 0;
}

void swap(struct Elem *i,struct Elem *j)
{
  i->prev->next=i->next;
  i->next->prev=i->prev;
  i->next=j->next;
  i->prev=j;
  j->next=i;
  i->next->prev=i;
}

void insertsort(struct Elem *element)
{
  struct Elem *start;
  struct Elem *current;
  start=element->next->next;
  while(start!=element){
    current=start->prev;
    while(current!=element && compare(start,current)==1) current=current->prev;
    swap(start,current);
    start=start->next;
  }
  start=start->next;

}
/*
void check_cond (void)
{
    int p = 0;
    if (current->next->v < current->next->next->v) {
        printf("%d ", current->next->v);
        goto exitcond;
    } else {
        if (current->next->v > current->next->next->v) {
            printf("%d ", current->next->next->v);
            p = 1;
            goto exitcond2;
        } else {
            p = 0;
            goto exitcond2;
        }
    }
    exitcond:;
    printf("%d", t1->next->next->v);
    return;
    exitcond2:;
    printf ("%d", t1->next->v);
    if (1 == p) return;
    printf(" %d", t1->next->next->v);
}*/

void mem_salvation (int marker)
{
    int arr[marker];
    int i;
    if (1 == marker) {
        free (current->next);
        free (current);
        return;
    }
    if (2 == marker) {
        free (t1->next->next);
        free (t1->next);
        free (t1);
        return;
    }
    if (2 < marker) {
        for (i = 0; i < marker; i++) arr[i] = t1->next;
        for (i = 0; i < marker; i++) free(arr[i]);
        free (t1);
    }
}


int main()
{
  int i, x;
  scanf("%d", &n);
  current = initlist();

  for(i = 0; i < n; i++){
        element = (struct Elem*)malloc(sizeof(struct Elem));
        scanf("%d", &x);
        element->v = x;
        Insert (current, element);
    }

    switch (n){
    case 0: return 0;
    case 1: {
        printf("%d",current->next->v);
        mem_salvation (n);
        return 0;
    }
    case 2: {
      if (current->next->v < current->next->next->v) {
            printf("%d ", current->next->v);
            printf("%d", current->next->next->v);
        } else {
            printf("%d ", current->next->next->v);
             printf("%d ", current->next->v);
        }

      free(current->next->next);
      free(current->next);
      free(current);

      return 0;
    }
    default: {
      insertsort(current);
      for(i = 0, print = current->next; i < n; i++, print = print->next) 
                printf("%d ",print->v);

      int f[n];

      for(i=0;i<n;i++) {
            f[i]=current;
            current=current->next;
      }

      for(i=0;i<n;i++) free(f[i]);

      free(current);
      return 0;
    }
  }
}