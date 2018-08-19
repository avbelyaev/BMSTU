.psect code gbl,rel,con,i,ro

;void add88( int88_x c, int88_x a, int88_x b )
;c = a + b
;typedef int int88_x[6];
;LE

madd::  mov R2, -(SP)
       mov R3, -(SP)
       mov R4, -(SP)
       mov R5, -(SP)

       mov 12(SP), R0  	;C
       mov 14(SP), R1  	;A
       mov 16(SP), R2  	;B

       mov #5, R3	;counter
       clr R4

loop:   mov R4, R5	; + prev carry
       clr R4		
       add (R1)+, R5	; + A[]
       adc R4
       add (R2)+, R5	; + B[]
       adc R4
       mov R5, (R0)+	;C[] = A[] + B[]
       sob R3, loop

;для хранения вектора используется нечетное число байт
;в последенем слове используется только байт с младшим адресом
;^^если такого условия нет, то опустить то, что ниже

       movb R4, R5	; + carry
       add (R1), R5    	; + A
       add (R2), R5    	; + B
       movb R5, (R0)

       mov 12(SP), R0
       mov (SP)+, R5
       mov (SP)+, R4
       mov (SP)+, R3
       mov (SP)+, R2

       return
.end


//Main.c (Decus C):

typedef int int88_x[6];

main()

{
		int88_x a, b, c;

		a[0] = 1; a[1] = 2; a[2] = 3; a[3] = 4; a[4] = 5; a[5] = 6;
		b[0] = 1; a[1] = 3; a[2] = 3; a[3] = 7; a[4] = 6; a[5] = 5;
		madd ( c, a, b );

		printf( "\r|%04x|%04x|%04x|%04x|%04x|%04x\r\n", 
			c[0], c[1], c[2], c[3], c[4], c[5] );

		return 0;
}

