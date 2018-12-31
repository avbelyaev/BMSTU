#!/bin/bash
flex lexer.l
bison -d parser.y
gcc -o test lex.yy.c parser.tab.c
