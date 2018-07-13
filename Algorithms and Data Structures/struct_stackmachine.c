#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define max(x, y) ((x)>(y) ? (x):(y))
#define min(x, y) ((x)<(y) ? (x):(y))

struct stack {
    int *data;
    int cap;
    int top;
};

int n;

char m_const[6] = "CONST";//push
char m_add[4] = "ADD";//+
char m_sub[4] = "SUB";//-
char m_mul[4] = "MUL";//*
char m_div[4] = "DIV";///
char m_max[4] = "MAX";
char m_min[4] = "MIN";
char m_neg[4] = "NEG";//sign-swap
char m_dup[4] = "DUP";//copy
char swag[5] = "SWAP";

void push (struct stack *s, int x)
{
    if (s->top == s->cap) printf("stack_overflow\n");
    else {
        s->data[s->top] = x;
        s->top++;
    }
}

int pop (struct stack *s)
{
    if (0 == s->top) {
        printf("devastation\n");
        return -1;
    }
    else {
        s->top--;
        int ret_val = s->data[s->top];
        return ret_val;
    }
}

void s_add (struct stack *s)
{
    int a = pop(s);
    int b = pop(s);
    push (s, a + b);
    //printf("add=%d\n", ret_val);
}

void s_sub (struct stack *s)
{
    int a = pop(s);
    int b = pop(s);
    push (s, a - b);
    //printf("sub=%d\n", ret_val);
}

void s_mul (struct stack *s)
{
    int a = pop(s);
    int b = pop(s);
    push (s, a * b);
    //printf("mul=%d\n", ret_val);
}

void s_div (struct stack *s)
{
    int a = pop(s);
    int b = pop(s);
    push (s, a / b);
    //printf("div=%d\n", ret_val);
}

void s_max (struct stack *s)
{
    int a = pop(s);
    int b = pop(s);
    push (s, max(a, b));
    //printf("max=%d\n", ret_val);
}

void s_min (struct stack *s)
{
    int a = pop(s);
    int b = pop(s);
    push (s, min(a, b));
    //printf("min=%d\n", ret_val);
}

void s_neg (struct stack *s)
{
    int a = pop(s);
    push (s, -a);
    //printf("neg=%d\n", ret_val);
}

void s_dup (struct stack *s)
{
    int a = pop(s);
    push (s, a);
    push (s, a);
    //printf("dup=%d\n", ret_val);
}

void s_swag (struct stack *s)
{
    int a = pop(s);
    int b = pop(s);
    push (s, a);
    push (s, b);
    //printf("swap\n");
}
void initstack (struct stack *s)
{
    s->data = (int*)malloc(n * sizeof(int));
    s->cap = n;
    s->top = 0;
}

int main()
{
    struct stack turn;
    int i, val;
    char s[6];

    scanf("%d", &n);

    initstack(&turn);

    for (i = 0; i < n; i++) {
        scanf("%s", s);
        if (0 == strcmp(s, m_const)) {
            scanf("%d", &val);
            push (&turn, val);
        }
        if (0 == strcmp(s, m_add)) {
            s_add (&turn);
        }
        if (0 == strcmp(s, m_sub)) {
            s_sub (&turn);
        }
        if (0 == strcmp(s, m_mul)) {
            s_mul (&turn);
        }
        if (0 == strcmp(s, m_div)) {
            s_div (&turn);
        }
        if (0 == strcmp(s, m_max)) {
            s_max (&turn);
        }
        if (0 == strcmp(s, m_min)) {
            s_min (&turn);
        }
        if (0 == strcmp(s, m_neg)) {
            s_neg (&turn);
        }
        if (0 == strcmp(s, m_dup)) {
            s_dup (&turn);
        }
        if (0 == strcmp(s, swag)) {
            s_swag (&turn);
        }
    }

    printf("%d", pop(&turn));

    free(turn.data);
    return 0;
}