#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAX_STRING_LEN 20
#define MAX_STRING_NUM 15

struct Elem {
    struct Elem *next;
    char *word;
};

char array[MAX_STRING_NUM][MAX_STRING_LEN], array2[10][10], s[35] = { 0 };
struct Elem *turn, *swp, *inl, *z, *bsa, *oper, *element;
int wnum, count;

int wcount (char *dest) {
    int i, n = 0;
    for (i = 0; i != strlen(dest); i++) if ((dest[i] != ' ') && (dest[i] != 9) && (dest[i+1] == ' ' || dest[i+1] == 13 || dest[i+1] == 9 || dest[i+1] == 0)) n++;
    wnum = n;
    return n;
}

struct Elem* InitSingleLinkedList ()
{
  struct Elem *first=(struct Elem*)malloc(sizeof(struct Elem));
  first->next=first;
  return first;
}


void insertafter(struct Elem *x, struct Elem *y)
{
  struct Elem *z;
  z=x->next;
  x->next=y;
  y->next=z;
}

void initarray (char *dest)
{
    int i = 0, j = 0, k = 0;
    while (i < strlen (dest)) {
        while (i < strlen(dest) && dest[i] == ' ') i++;
        while (i < strlen(dest) && dest[i] != ' ') {
            array[j][k] = dest[i];
            k++;
            i++;
        }
        j++;
        k = 0;
    }
}

int compare(struct Elem *a, struct Elem *b)
{
  if (strlen(a->word)>strlen(b->word)) return 1;
  return 0;
}

void swap(struct Elem *a, struct Elem *b)
{
  struct Elem *s;

  s=(struct Elem*)malloc(sizeof(struct Elem));
  s->word=a->word;
  a->word=b->word;
  b->word=s->word;
  free(s);
}

void result (struct Elem *list)
{
    int i;
    bsa = list;
    for(i = 0; i < wnum; i++) {
        printf("%s ", bsa->word);
        bsa = bsa->next;
    }
}

struct Elem *bsort(struct Elem *list)
{
  int i, j;
  for(i = wnum - 1; i > 0; i--)
          for(j = 0, bsa = list; j < i; j++, bsa = bsa->next)
          	if (1 == compare(bsa, bsa->next)) swap (bsa, bsa->next);
  		bsa = list;
  	for(i = 0; i < wnum; i++) {
        printf("%s ", bsa->word);
        bsa = bsa->next;
    }
}

int main()
{
  int i, j;
  char s[1000];
  gets(s);
  wnum = wcount(s);
  initarray(s);

  struct Elem *current;
  current=InitSingleLinkedList();

  struct Elem *element;

  for(i = wnum; i >= 0; i--){
    element = (struct Elem*)malloc(sizeof(struct Elem));
    insertafter (current, element);
    element->word = array[i];
  }

  bsort(current->next);
  //mem_salvation();
  int n = wnum;
  current = current->next;
  int freed [n+1];

  for(i = 0; i <= n; i++) {
        freed[i] = current;
        current = current->next;
  }

  for(i = 0; i <= n; i++) free (freed[i]);

  free(current);
  return 0;
}