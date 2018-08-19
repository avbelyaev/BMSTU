.psect code,gbl,i,con,ro

;lsh72(b,a,n)
;b = a << n
;typedef int int72_x[5]
;LE

lsh72:: mov R5, -(SP)
        mov SP, R5
        mov R4, -(SP)
        mov R3, -(SP)
        mov R2, -(SP)           ;R2,R3,R4,R5,ret...
 
	mov 4(R5), R1		;B
        mov 6(R5), R3		;A
        
;copy A to B
;<-- better add here cycle -->
        mov R1, R4
        mov (R3)+, (R4)+
        mov (R3)+, (R4)+
        mov (R3)+, (R4)+
        mov (R3)+, (R4)+
        movb (R3)+, (R4)+

	mov 10(R5), R0		;N
        beq EXIT

loop:   mov R1, R4        
        clr R2
        clr R5
        asl (R4)+       ;[0]<<
        adc R2          ;R2=carry[0]

        asl (R4)        ;[1]<<
        adc R5          ;R5=carry1
        add R2, (R4)+
        clr R2          

        asl (R4)        ;[2]<<
        adc R2          ;R2=carry2
        add R5, (R4)+
        clr R5         

        asl (R4)        ;[3]<<
        adc R5          ;R5=carry3
        add R2, (R4)+
        clr R2          ;R2=0

        movb (R4), R2   ;R2=b[last]
        asl R2          ;b[last]<<
        movb R2, (R4)   ;R4=b[last]<<
        add R5, (R4)    ;+carry3

        sob R0, loop

EXIT:   mov (SP)+, R2
        mov (SP)+, R3
        mov (SP)+, R4
        mov (SP)+, R5

        return
.end

