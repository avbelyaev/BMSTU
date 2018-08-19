#include <iostream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <algorithm>

using namespace std;

int var, tag, num, indx;
string image, number;

struct stack {
    int data[100];
    int top;
};

void InitStack (struct stack *s)
{
    s->top = 0;
}
int pop (struct stack *s)
{
    return s->data[--s->top];
}
void push (struct stack *s, int x)
{
    s->data[s->top++] = x;
}
int show (struct stack *s)
{
    int retval = s->data[--s->top];
    s->data[s->top++] = retval;
    return retval;
}

void parse_E (string str, struct stack *s);
void parse_E1 (string str, struct stack *s);
void parse_T (string str, struct stack *s);
void parse_T1 (string str, struct stack *s);
void parse_F (string str, struct stack *s);

int prior (const char *str, int x)
{
    if (str[x] == '+' || str[x] == '-') return 1;
    if (str[x] == '*' || str[x] == '/') return 2;
}

int sign (const char *str, int x)
{
    if (str[x] == '+') return 1;
    if (str[x] == '-') return 2;
    if (str[x] == '*') return 3;
    if (str[x] == '/') return 4;
}

int myatoi(string string_str)
{
    const char * char_str = string_str.c_str();
    int val = 0;
    while( *char_str ) val = val*10 + (*char_str++ - '0');
    return val;
}

void NextLexem (string str)
{
    int jndx = 0;
    while (indx < str.size()) {
        if ('_' != str[indx]) {
            if (str[indx] == '+') {
                image = '+';
            } else if (str[indx] == '-') {
                image = '-';
                } else if (str[indx] == '*') {
                    image = '*';
                    } else if (str[indx] == '/') {
                        image = '/';
                        } else if (str[indx] == '(') {
                            image = '(';
                            } else if (str[indx] == ')') {
                                image = ')';
                            }
            tag = 1;
            if (str[indx] == 'x') {
                tag = 2;
                image = 'x';
            }
            jndx = indx;
            if ((str[indx] >= '0') && (str[indx] <= '9')) {
                tag = 3;
                number = "";
                while ((jndx < str.size()) && (str[jndx] >= '0') && (str[jndx] <= '9')) {
                    number += str[jndx];
                    jndx++;
                }
                indx = jndx - 1;
            }
        } else {
            indx++;
            continue;
        }
        indx++;
        break;
    }

}

void parse_E (string str, struct stack *s)
{
    parse_T (str, s);
    parse_E1 (str, s);
}

void parse_T (string str, struct stack *s)
{
    parse_F (str, s);
    parse_T1 (str, s);
}

void parse_E1 (string str, struct stack *s)
{
    if ((tag == 1) && (image == "-")) {
        NextLexem(str); //next lexem for <T>
        parse_T(str, s);   //<T>
        int arg1 = pop(s), arg2 = pop(s);
        push(s, arg2 - arg1);
        parse_E1(str, s);  //<E1>
    } else if ((tag == 1) && (image == "+")) {
        NextLexem(str);
        parse_T(str, s);
        push(s, pop(s) + pop(s));
        parse_E1(str, s);
        }
}

void parse_T1 (string str, struct stack *s)
{
    if ((tag == 1) && (image == "/")) {
        NextLexem(str); //next lexem for <F>
        parse_F(str, s);   //<F>
        int arg1 = pop(s), arg2 = pop(s);
        push(s, arg2 / arg1);
        parse_T1(str, s);  //<T1>
    } else if (tag == 1 && image == "*") {
        NextLexem(str);
        parse_F(str, s);
        push(s, pop(s) * pop(s));
        parse_T1(str, s);
        }
}

void parse_F (string str, struct stack *s)
{
    if (tag == 1) {
        if (image == "-") {
            NextLexem(str);
            parse_F(str, s);
            push(s, pop(s) * (-1));
        } else if (image == "(") {
            NextLexem(str); //next lexem for <E>
            parse_E(str, s);   //( <E> )
            NextLexem(str);
        }
    } else if (tag == 2) {
            push(s, var);
            NextLexem(str);
            } else if (tag == 3) {
                    push(s, myatoi(number));
                    NextLexem(str);
            }
}

string lucky (string expr)
{
    const char * str = expr.c_str();
    int i, delnum = 0;
    for (i = 0; i <= strlen(str); i++)
        if (str[i] == '(' || str[i] == ')') expr.erase(i - delnum++, 1);
    return expr;
}

string del_unary_n_extrenal (string expr)
{
    const char * str = expr.c_str();
    int i, j, buf = 0;
    for (i = strlen(str); i >= 0; i--) {
        if ((str[i] == '-' || str[i] == '+' || str[i] == '*' || str[i] == '/') && str[i+1] == '(' && str[i+3] == ')') {
            expr.erase(i+3, 1);
            expr.erase(i+1, 1);
        }
    }

    const char * str_no_unary = expr.c_str();
    if (str_no_unary[0] == '(' && str_no_unary[strlen(str_no_unary)-1] == ')') {
        for (i = 0; i <= strlen(str_no_unary); i++) {
            if (str_no_unary[i] == '(') buf++;
            if (str_no_unary[i] == ')') buf--;
            if (buf == 0) break;
            if (i == strlen(str_no_unary) - 2 && buf == 1) {
                expr.erase(strlen(str_no_unary)-1,1);
                expr.erase(0, 1);
            }
        }
    }
    return expr;
}

void del_brackets (string expr)
{
    const char * str = expr.c_str();
    int i, j, k, buf, prior1, prior2, sign1, sign2, delnum = 0;
    for (i = 0; i < expr.size(); i++) {
        if (str[i] == '+' || str[i] == '-' || str[i] == '*' || str[i] == '/') {
            prior1 = prior(str, i);
            sign1 = sign(str, i);
            if (str[i-1] == ')') {
                    buf = 0;
                    j = i - 1;
                    while (j >=0) {
                        if (str[j] == '(') buf++;
                        if (str[j] == ')') buf--;
                        if (buf == 0) break;
                        j--;
                    }
                    if (buf == 0) {
                            for (k = i - 1; k >= j; k--) {
                                    prior2 = 2;
                                if (str[k] == '+') {
                                    prior2 = 1;
                                    sign2 = 1;
                                    break;
                                }
                                if (str[k] == '-') {
                                    prior2 = 1;
                                    sign2 = 2;
                                    break;
                                }
                                if (str[k] == '*' || str[k] == '/') prior2 = 2;
                            }
                            if (prior2 == 2 || (prior1 == 1 && sign2 == 1) || (prior1 == 1 && sign2 == 2)) {
                                    expr.erase(j, 1);
                                    delnum++;
                                    expr.erase(i-1-delnum++,1);
                            }
                    }
                }
            if (str[i+1] == '(') {
                    buf = 0;
                    j = i + 1;
                    while (j <= strlen(str)) {
                        if (str[j] == '(') buf++;
                        if (str[j] == ')') buf--;
                        if (buf == 0) break;
                        j++;
                    }
                    if (buf == 0) {
                            for (k = i+1; k <= j; k++) {
                                    prior2 = 2;
                                if (str[k] == '+') {
                                    prior2 = 1;
                                    sign2 = 1;
                                    break;
                                }
                                if (str[k] == '-') {
                                    prior2 = 1;
                                    sign2 = 2;
                                    break;
                                }
                                if (str[k] == '*' || str[k] == '/') prior2 = 2;
                            }
                            if (sign1 == 1 || (sign1 == 2 && prior2 == 2) || (prior1 == 2 && prior2 == 2)) {
                                    expr.erase(i+1-delnum++, 1);
                                    expr.erase(j-delnum++,1);
                            }
                    }
                }
        }
    }
    cout << expr << endl;
}

void parse (string expr1)
{
    struct stack sm1;
    indx = 0;
    number = "";
    InitStack(&sm1);
    NextLexem(expr1);
    parse_E(expr1, &sm1);
    cout << show(&sm1) << endl;

    struct stack sm2;
    indx = 0;
    number = "";
    InitStack(&sm2);
    string expr2 = lucky(expr1);
    NextLexem(expr2);
    parse_E(expr2, &sm2);

    if (show(&sm1) == show(&sm2)) {
        cout << expr2 << endl;
        return;
    }
    del_brackets (del_unary_n_extrenal(expr1));
}

int main()
{
    string expr;
    cin >> expr;
    cin >> var;
    parse(expr);
    return 0;
}

