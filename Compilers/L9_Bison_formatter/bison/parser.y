
%{
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lexer.h"

#define TYPENAME_FLOAT "float"
#define TYPENAME_INT "int"

#define MEM(size) ((char *)malloc(sizeof(char) * (size + 1)))
%}

%define api.pure
%locations
%lex-param {yyscan_t scanner}
%parse-param {yyscan_t scanner}

%union {
	char *string;
	struct number {
		char *typeName;
		union numData {
			int numInt;
			float numFloat;
		} numData;
	} number;
}

%token <string> String
%token <number> Number
%token ValTrue ValFalse ValNull
%left ObjBeg ObjEnd ArrBeg ArrEnd SymLBreak SymTab
%left SymComma
%left SymColon
%start Start

%type<string> Object Array ObjBeg ObjEnd ArrBeg ArrEnd Members Pair Elements Value x

%{
int yylex(YYSTYPE *yylval_param, YYLTYPE *yylloc_param , yyscan_t scanner);
void yyerror(YYLTYPE *yylloc, yyscan_t scanner, char *msg);
%}

//	$$			$1		  $2			$3			   $4		 $5		
//  Members		ObjBeg 	  x 			Members 	   ObjEnd 	 x 
//				"{\t"	  %s 			%s 			   " }"		 %s
//  			3+		  strlen($2)+	strlen($3)+	   2+		 strlen($5)

%%

Start:		Object 						{ printf("%s\n",$1); }
	;

Object:		ObjBeg x ObjEnd x			{ $$ = MEM(1+strlen($2)+1+strlen($4)); sprintf($$, "{%s}%s", $2,$4); }
	|		ObjBeg x Members ObjEnd x	{ $$ = MEM(3+strlen($2)+strlen($3)+2+strlen($5)); sprintf($$,"{\t%s%s }%s", $2,$3,$5); }
	;

Members:	Pair 						{ $$ = $1; }
	|		Pair SymComma x Members		{ $$ = MEM(2+strlen($1)+2+strlen($3)+strlen($4)); sprintf($$,"\t%s, %s%s", $1,$3,$4); }
	;

Pair:		String x SymColon x Value	{ $$ = MEM(strlen($1)+strlen($2)+3+strlen($4)+strlen($5)); sprintf($$,"%s%s : %s%s", $1,$2,$4,$5); }
	;

Array:		ArrBeg x ArrEnd	x			{ $$ = MEM(1+strlen($2)+1+strlen($4)); sprintf($$,"[%s]%s", $2,$4); }
	|		ArrBeg x Elements ArrEnd x	{ $$ = MEM(2+strlen($2)+strlen($3)+2+strlen($5)); sprintf($$,"[ %s%s ]%s", $2,$3,$5); }
	;

Elements:	Value						{ $$ = $1; }
	|		Value SymComma x Elements	{ $$ = MEM(strlen($1)+2+strlen($3)+strlen($4)); sprintf($$,"%s, %s%s",$1,$3,$4); }
	;

Value:		Object						{ $$ = $1; }
	|		Array						{ $$ = $1; }
	|		String x					{ $$ = MEM(strlen(yylval.string)+strlen($2)); sprintf($$,"%s%s", yylval.string, $2); }
	|		Number x					{ 
											if (0 == strcmp(yylval.number.typeName, TYPENAME_INT)) {
												$$ = malloc(sizeof(int) * 1 + sizeof(char*) * (strlen($2)+1)); 
												sprintf($$, "%d%s", yylval.number.numData.numInt, $2);
											} else {
												$$ = malloc(sizeof(float) * 1 + sizeof(char*) * (strlen($2)+1)); 
												sprintf($$, "%f%s", yylval.number.numData.numFloat, $2);
											}
		 								}
	|		ValTrue	x					{ $$ = MEM(4+strlen($2)); sprintf($$, "true%s",  $2); }
	|		ValFalse x					{ $$ = MEM(5+strlen($2)); sprintf($$, "false%s", $2); }
	|		ValNull	x					{ $$ = MEM(4+strlen($2)); sprintf($$, "null%s",  $2); }
	;

x:										{ $$ = ""; }
	|		SymLBreak x					{ $$ = MEM(2+strlen($2)); sprintf($$,"\n%s", $2); } 
	|		SymTab x					{ $$ = MEM(2+strlen($2)); sprintf($$,"\t%s", $2); } 
	;

%%    

int main()
{
	char * buffer = 0;
	long length;
	FILE * f = fopen ("input.txt", "rb");

	if (f) {
		fseek (f, 0, SEEK_END);
		length = ftell (f);
		fseek (f, 0, SEEK_SET);
		buffer = malloc (length);
		if (buffer) 
			fread (buffer, 1, length, f);
		fclose (f);
	}

	if (buffer) {
		printf("===============================================\n");
		printf("----------------- INPUT JSON ------------------\n");
		printf("===============================================\n");
		printf("%s\n", buffer);
		printf("===============================================\n");
		printf("---------------- OUTPUT JSON ------------------\n");
		printf("===============================================\n");

		//Weak formatter
		yyscan_t scanner;
		struct Extra extra;
		init_scanner(buffer, &scanner, &extra);
		yyparse(scanner);
		destroy_scanner(scanner);

		free(buffer);
	}
	return 0;
}
