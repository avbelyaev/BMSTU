public interface StackMachine 
{ 
        // Положить константу в стек 
        void pushConst(int value); 
 
        // Положить в стек значение переменной x 
        void pushVar(); 
 
        // Изменить знак числа на вершине стека 
        void neg(); 
 
        // Бинарные арифметические операции. 
        // Каждая операция снимает два операнда со стека 
        // и кладёт на стек результат 
        void add();  // сложение 
        void sub();  // вычитание 
        void mul();  // умножение 
        void div();  // деление 
}

//========================================================================================

    public class ParsingDirector {
            private StackMachine stack;
            private String number;
            private char image;
            private int tag, i, num;
    /*LL1 - грамматика:
    <E>  ::= <T> <E’>.
    <E’> ::= + <T> <E’> | - <T> <E’> | .
    <T>  ::= <F> <T’>.  
    <T’> ::= * <F> <T’> | / <F> <T’> | .
    <F>  ::= <number> | <var> | ( <E> ) | - <F>.
    */    
            public ParsingDirector (StackMachine stack) { this.stack = stack; }
           
            public void parse (String s) {
                    i = 0;
                    NextLexem(s);
                    parse_E(s);     //<E>  ::= <T> <E’>.
            }
           
            private void NextLexem (String str) {
                    int j = 0;
                    while (i < str.length()) {
                            if (str.charAt(i) != ' ') {
                                    if (str.charAt(i) == '+') {
                                            image = '+';
                                    } else if (str.charAt(i) == '-') {
                                            image = '-';
                                    } else if (str.charAt(i) == '*') {
                                                    image = '*';
                                            } else if (str.charAt(i) == '/') {
                                                            image = '/';
                                                    } else if (str.charAt(i) == '(') {
                                                                    image = '(';
                                                            } else if (str.charAt(i) == ')') {
                                                                            image = ')';
                                                            }
                            tag = 1;
                            if (str.charAt(i) == 'x') {
                                    tag = 2;
                                    image = 'x';
                            }
                            j = i;
                            if ((str.charAt(i) >= '0') && (str.charAt(i) <= '9')) {
                                    tag = 3;
                                    number = "";
                                    while ((j < str.length()) && (str.charAt(j) >= '0') && (str.charAt(j) <= '9')) {
                                            number += str.charAt(j);
                                            j++;    //move in number
                                    }
                                    i = j-1;
                            }
                            } else {
                                    i++;
                                    continue;
                            }
                            i++;
                            break;
                    }
            }
           
            //<E>  ::= <T> <E’>.
            private void parse_E (String str) {
                    parse_T (str);
                    parse_E1 (str);
            }
           
            //<T>  ::= <F> <T’>.
            private void parse_T (String str) {
                    parse_F (str);
                    parse_T1 (str);
            }
           
            //<E’> ::= + <T> <E’> | - <T> <E’> | .
            private void parse_E1 (String str) {
                    if (tag == 1 && image == '-') {
                            NextLexem(str); //next lexem for <T>
                            parse_T(str);   //<T>
                            stack.sub();    //- <T> <E1>
                            parse_E1(str);  //<E1>
                    } else if (tag == 1 && image == '+') {
                            NextLexem(str);
                            parse_T(str);
                            stack.add();
                            parse_E1(str);
                    }
            }
           
            //<T’> ::= * <F> <T’> | / <F> <T’> | .
            private void parse_T1 (String str) {
                    if (tag == 1 && image == '/') {
                            NextLexem(str); //next lexem for <F>
                            parse_F(str);   //<F>
                            stack.div();    // / <F> <T1>
                            parse_T1(str);  //<T1>
                    } else if (tag == 1 && image == '*') {
                            NextLexem(str);
                            parse_F(str);
                            stack.mul();
                            parse_T1(str);
                    }
            }
           
            //<F>  ::= <number> | <var> | ( <E> ) | - <F>.
            private void parse_F (String str) {
                    if (tag == 1) {
                            if (image == '-') {
                                    NextLexem(str);
                                    parse_F(str);
                                    stack.neg(); //- <F>
                            } else if (image == '(') {
                                    NextLexem(str); //next lexem for <E>
                                    parse_E(str);   //( <E> )
                                    NextLexem(str);
                            }
                    } else {
                            if (tag == 2) {
                                    stack.pushVar();        //<var>
                                    NextLexem(str);
                            } else {
                                    if (tag == 3) {
                                            num = Integer.parseInt(number);
                                            stack.pushConst(num);   //<number>
                                            NextLexem(str);
                                    }
                            }
                    }
            }
    }

//========================================================================================

=> 5 
=> x*(x+2)
<= 35

=> 10 
=> 2 * -x + 100
<= 80

