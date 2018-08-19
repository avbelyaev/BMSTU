.586
_TEXT   segment byte public 'CODE' use16
assume cs:_TEXT, ds:nothing

org     100h	
			
start:
        cli				;Clear Interrupt Flag
        lss     SP, dword ptr STKPTR
        sti				;Set Interrupt Flag

	call  	init
	call 	stop
  
STKPTR  dw	0FFFEh, 09000h

init:	
	mov  dx, 0C000h  		;начинаем с видеокарты
cycle: 	mov  DS, dx 			
  	mov  ax, 80h			;в адр пространстве модуль занимает диапазон, кратный 2 кб
  	cmp  word ptr DS:[0], 0AA55h	;проверка нначала модуля
  	jnz  nxt			
  
	call scanbios

	;xor ax, ax
  	movzx  ax, byte ptr DS:[2]	;размер модуля (3й байт DS) -> ax

;округление вверх
  	add  al, 3h			;число блоков по 512 должно быть кратно 2048/512 (если в конце числа 00 - кратно 4)
  	and  al, 0FCh			;накладываем маску (11111100)

;переход физ -> лог адр
	shl  ax, 5  			;N_блоков * 512 / 16 эквивалентно сдвигу на 5
					;512 = 2^9 - размер блока
					;0x10 = 2^4 - размер параграфа

nxt: 	add  dx, ax			;+2кб
  	cmp  dx, 0F000h			;проверка границы области модулей
  	jb  cycle			
	ret
  
scanbios proc near
  	cld				;clear direction flag
  	xor     si, si
 	xor     cx, cx
  	mov     ch, DS:[2]		
  	xor     bl, bl

chcksm: lodsw				;переписывает слово, чей адрес в памяти определяется парой регистров DS: (E)SI, в регистр ax
        add     al, ah			
        add     bl, al			
        dec     cx
        jnz     short chcksm

        or      bl, bl			;сумма по модулю 256 всех байт от xxxxx+0 до xxxxx+size*512 должна быть равна 0
        jnz     short skip		

        pusha
        push    ds
        push    es
        push    fs
        push    gs

        push    CS
        push    offset __ret
        push    DS
        push    3h
        retf
__ret:
        pop     gs
        pop     fs
        pop     es
        pop     ds
        popa
skip:	ret     0
scanbios endp

stop	proc	near
	cli
	hlt
	jmp	short stop
stop	endp

org	0FFF0h
	db	0EAh		; JMP FAR
	dw	offset start	; offset
	dw	0F000h		; segment


org	0FFFEh
	dw	99FCh		; PC 

_TEXT	ends
end	start

