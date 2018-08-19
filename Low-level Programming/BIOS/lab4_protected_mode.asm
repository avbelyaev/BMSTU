.586p

include mymacro.asm

;========================================================================
;-----------------------ADDRESSES-AND-VARIABLES-------------------------
;========================================================================

;-----------------------------DESCRIPTORS-------------------------------
DCODE		equ 8			;дескриптор кода
DDATA		equ 10h			;дескриптор данных
DSTACK		equ 18h			;дескриптор стека
DIOAPIC		equ 20h
DLAPIC		equ 28h

STARTGDT	equ 50h			;начальный адрес GDT
STARTIDT	equ 60h			;начальный адрес IDT

VIDEOB		equ	6000h		;база видеопамяти
STACK		equ	9F000h		;база стека
CODE		equ	0F0000h		;сегмент кода
;--------------------------------UART-----------------------------------
COM1		equ 	3F8h
COM2		equ 	2F8h
COM3		equ 	3E8h
COM4		equ 	2E8h

PORT		equ	000h
THR		equ 	000h		;write only
RBR		equ 	000h		;read only
IER		equ	001h
DLL		equ	000h
DLH		equ	001h
IIR		equ 	002h		;read only
FCR		equ	002h		;write only
LCR		equ	003h
MCR		equ	004h
LSR		equ	005h		;read only
MSR		equ	006h		;read only
SR		equ 	007h		;scratch
;-------------------------------BUFFER-----------------------------------
SBUFADR		equ	0h		;адрес текущего начала буфера
EBUFADR		equ	2h		;адрес текущего конца буфера
SBUF		equ 	6h		;начало буфера
EBUF		equ	3Ah		;конец буфера
;-----------------------------CONTROLLERS--------------------------------
TESTOK 		equ	0AAh		
ACK			equ	0FAh
;---------------------------IOAPICS-LAPICS-------------------------------
IOAPICB		equ 0FEC00000h	;адрес иоапика
LAPICB		equ 0FEE00000h	;адрес лапика

;========================================================================
;--------------------------------START----------------------------------
;========================================================================
_TEXT32 segment byte public 'code' use32
assume cs:_TEXT32, ds:nothing
_TEXT32START = $
		org		100h
start32:
		mov ax,	DDATA					;перезагружаем сегментные регистры.
		mov ds, ax
		mov es, ax
		mov fs, ax
		mov gs, ax
		mov ax, DSTACK
		mov ss, ax
		mov esp, 0FFCh

		call	permit_NMI				;NMI - Not Masked Interruptions
		sti
		call	init_PIC
		call	init_IOAPIC
		call	init_LAPIC
		call 	init_UART
		call	FINAL_CYCLE
		call 	stop

;------------------------------------------------------------------------
permit_NMI:
		in al, 70h
		and al, 7Fh
		out 70h, al
		ret


;========================================================================
;-------------------------------TABLES----------------------------------
;========================================================================
;--------------------------- KEYBOARD INITIALIZATION TABLE -------------
KBD_INIT	db	0AAh,	064h,	1,		055h,				
			0ABh,	064h,	1,		000h,				
			0AEh,	064h,	0,						
			0FFh,	060h,	2,	ACK, TESTOK				
		db	0EEh			;stop-symbol				
;------------------------------------------------------------------------

;----------------------------- COM INTERVIEW TABLE ---------------------    
COM_TABLE	dw	COM1,	COM2,	COM3,	COM4
;------------------------------------------------------------------------

;----------------------------- UART INITIALIZATION TABLE ---------------
COM_INIT	db	10000000b,	LCR,	
			60h,		DLL,	
			0,		DLH,	
			00000010b,	LCR,	
			0,		IER,
			00000110b,	FCR,		
			00001000b,	MCR									
		db	0FEh			;stop-symbol
;1200 baud, 7 bit, 1 stop bit, no parity
;------------------------------------------------------------------------

move_out macro	port, value
	mov al, value
	out port, al
endm

mov_segm macro	sreg, value, treg:vararg
	_treg equ ax
	for t, <treg>
		t
	endm
	mov _treg, value
	mov sreg, _treg
endm

req_resp macro port, rq, rp:vararg
	mov dx, port
	mov bl, rq
	call request
	for rsp, <rp>
		mov bl, rsp
		call respond
		cmp al, bl
		;jnz stop
	endm
endm

send_com macro	port, value
		push dx
		mov		dx, port
		mov		al, value
		out		dx, al
		pop dx
endm
;========================================================================
;---------------------------------init_PIC------------------------------
;========================================================================
init_PIC proc near
;начальный сброс обоих контроллеров

	mov al, 20h		;1
	out 20h, al
	mov al, 20h		;2
	out 0A0h, al

	mov al, 11h		;1
	out 20h, al
	mov al, 11h		;2
	out 0A0h, al

	mov al, 08h		;1
	out 21h, al
	mov al, 70h		;2
	out 0A1h, al

	mov al, 04h		;1
	out 21h, al
	mov al, 02h		;2
	out 0A1h, al
	
	mov al, 05h		;1
	out 21h, al
	mov al, 01h		;2
	out 0A1h, al

;разршение выбранных IRQ
;PIC0
	mov al, 0FBh 	;1111 1011
	out 21h, al
;PIC1
	mov al, 0FFh	
	out 0A1h, al
	ret
init_PIC endp



;========================================================================
;---------------------------------COUTS---------------------------------
;========================================================================
space proc near
		mov al, 20h
		call pushbuf
		ret
space endp

cout_kb proc near
		mov al, 6Bh				;"kb:"
		call pushbuf
		mov al, 62h
		call pushbuf
		mov al, 3Ah
		call pushbuf
		ret
cout_kb endp

cout_mf proc near
		mov al,	'M'				;"MF"
		call pushbuf
		mov al, 'F'
		call pushbuf
		ret
cout_mf endp

cout_at proc near
		mov al, 'A'				;"AT"
		call pushbuf
		mov al, 'T'
		call pushbuf
		ret
cout_at endp

cout_xt proc near
		mov al, 'X'				;"XT"
		call pushbuf
		mov al, 'T'
		call pushbuf
		ret
cout_xt endp

cout_minus proc near			;"-"
		mov al, '-'
		call pushbuf
		ret
cout_minus endp

cout_id proc near
		mov al, 69h				;"id:"
		call pushbuf
		mov al, 64h
		call pushbuf
		mov al, 3Ah
		call pushbuf
		ret
cout_id endp

cout_ms proc near				;"  ms:"	
		call space
		call space
		mov al, 'm'
		call pushbuf
		mov al, 73h
		call pushbuf
		mov al, 3Ah
		call pushbuf
		ret
cout_ms endp

cout_ps2 proc near
		mov al, 'p'
		call pushbuf
		mov al, 's'
		call pushbuf
		mov al, '/'
		call pushbuf
		mov al, '2'
		call pushbuf
		ret
cout_ps2 endp

cout_im1 proc near
		mov al, 'i'
		call pushbuf
		mov al, 'm'
		call pushbuf
		mov al, '1'
		call pushbuf
		ret
cout_im1 endp

cout_im2 proc near
		mov al, 'i'
		call pushbuf
		mov al, 'm'
		call pushbuf
		mov al, '2'
		call pushbuf
		ret
cout_im2 endp

cout_uart proc near				;"  uart:"
		call space
		call space
		mov al, 'u'					
		call pushbuf
		mov al, 'a'					
		call pushbuf
		mov al, 'r'					
		call pushbuf
		mov al, 't'					
		call pushbuf
		mov al, 3Ah
		call pushbuf
		ret
cout_uart endp

cout_16 proc near				;"16"
		mov al, '1'
		call pushbuf
		mov al, '6'
		call pushbuf
		ret
cout_16 endp

cout_bps proc near
		mov al, 'b'
		call pushbuf
		mov al, '/'
		call pushbuf
		mov al, 's'
		call pushbuf
		ret
cout_bps endp

cout_00 proc near
		mov al, '0'
		call pushbuf
		mov al, '0'
		call pushbuf
		ret
cout_00 endp

cout_type16750 proc near		;"16750 "
		call cout_16
		mov al, '7'
		call pushbuf
		mov al, '5'
		call pushbuf
		mov al, '0'
		call pushbuf
		call space
		ret
cout_type16750 endp

cout_type16550A proc near		;"16550A "
		call cout_16
		mov al, '5'
		call pushbuf
		mov al, '5'
		call pushbuf
		mov al, '0'
		call pushbuf
		mov al, 'A'
		call pushbuf
		call space
		ret
cout_type16550A endp

cout_type16550 proc near		;"16550 "
		call cout_16
		mov al, '5'
		call pushbuf
		mov al, '5'
		call pushbuf
		mov al, '0'
		call pushbuf
		call space
		ret
cout_type16550 endp

cout_type16450 proc near		;"16450 "
		call cout_16
		mov al, '4'
		call pushbuf
		mov al, '5'
		call pushbuf
		mov al, '0'
		call pushbuf
		call space
		ret
cout_type16450 endp

cout_type8250 proc near			;"8250 "
		mov al, '8'
		call pushbuf
		mov al, '2'
		call pushbuf
		mov al, '5'
		call pushbuf
		mov al, '0'
		call pushbuf
		call space
		ret
cout_type8250 endp

cout_2400bps proc near			;"2400b/s "
		mov al, '2'
		call pushbuf
		mov al, '4'
		call pushbuf
		call cout_00
		call cout_bps
		call space
		ret
cout_2400bps endp

cout_9600bps proc near			;"9600b/s "
		mov al, '9'
		call pushbuf
		mov al, '6'
		call pushbuf
		call cout_00
		call cout_bps
		call space
		ret
cout_9600bps endp

cout_19200bps proc near			;"19200b/s "
		mov al, '1'
		call pushbuf
		mov al, '9'
		call pushbuf
		mov al, '2'
		call pushbuf
		call cout_00
		call cout_bps
		call space
		ret
cout_19200bps endp

cout_115200bps proc near			;"115200b/s "
		mov al, '1'
		call pushbuf
		mov al, '1'
		call pushbuf
		mov al, '5'
		call pushbuf
		mov al, '2'
		call pushbuf
		call cout_00
		call cout_bps
		call space
		ret
cout_115200bps endp
;========================================================================
;---------------------------------init_IOAPIC---------------------------
;========================================================================
;1.IOAPIC обычно отображается на адрес IOAPICBASE 0xFEC00000 (адрес может быть несколько изменен средствами северного моста)
;2.IOZPIC представлен двумя регистрами: "селектор" (selector) по адресу IOAPICBASE и "окно" (window) по адресу IOAPICBASE+0x10
;3.При работе с IOAPIC в селектор записывается номер регистра IOAPIC, а из окна производится чтение или запись регистра
;4.Термин "Reserved", определяющий те или иные биты регистров IOAPIC обозначает не "MBZ", а то, что содержимое этих разрядов нельзя изменять, 
;при необходимости изменить регистр надо сначала прочитать прежнее значение, потом только нужные биты обнулить или установить, а затем записать обратно.
.586p
init_IOAPIC proc near
	push	DS
	push	cx
	smov	DS, DIOAPIC

;0 линия, 10 рег
	mov		dword ptr DS:[0h], 10h			;номер регистра отвечающего за 0-ю линию -> в адрес селектора (0-я линия - 10,11 рег, 23-я линия - 3E,3F рег)
	mov 	eax, dword ptr DS:[10h]			;работа с регистром в окне
  	and 	eax, 0FFFE60FFh							;Физическая адресация LAPIC, обнуляем Trigger Mode(15), Pin Polarity(13), возбуждение прерывания передним фронтом (состояние 1) (т.к. настройка irq0)
  	or 		eax, 0700h								;Delivery Mod = 111 = ExtInt - Causes the processor to respond the intr as if the intr originated in an externally connected intr controller								
  	mov 	dword ptr DS:[10h], eax	
;0 линия, 11 рег
  	mov 	dword ptr DS:[0h], 011h	
  	mov 	eax, dword ptr DS:[10h]	
  	and 	eax, 0FFFFFFh							;Обнулить целевой ID лапика приёмника (0-й lapic единственного процессора, который эмулирует бочс)
  	mov 	dword ptr DS:[10h], eax	

  	mov 	ecx, 012h								
loop_masking:
  	mov 	DS:[0h], ecx					;номер очередного регистра в селектор
  	mov 	eax, dword ptr DS:[10h]	
  	or 		eax, 010000h							;16-й бит MASK
  	mov 	dword ptr DS:[10h], eax	
  	add 	ecx, 2									
  	cmp 	ecx, 040h								
  	jnz 	loop_masking

  	pop 	cx								
  	pop 	DS					
  	ret
init_IOAPIC endp

;========================================================================
;---------------------------------init_LAPIC----------------------------
;========================================================================
init_LAPIC proc near
		push DS
		push cx
		mov	ecx, 01Bh
		rdmsr
		or	eax, 0800h
		wrmsr
		smov DS, DLAPIC

		;mov DS:[040h], offset spur_intr				;обработка фиктивных прерываний
		;mov DS:[042h], cs
		mov dword ptr DS:[0F0h], 0110h			;apicen = ApicEnabled и бит разрешения -> в дескриптор spurious вектора 
		or dword ptr DS:[350h], 010000h		;максирование LINT0 
		or dword ptr DS:[360h], 010000h		;маскирование LINT1
		;mov DS:[080h], offset interrupt				
		;mov DS:[082h], cs
		
		mov	eax, dword ptr DS:[320h]			;LVT[0] - Timer
		and	eax, 06F00h
		or	eax, 020020h							
		mov dword ptr DS:[320h], eax
		or dword ptr DS:[3E0h], 0Bh			;Timer Divide Config: 1011, d0=1, d1=1, d3=1 => делитель частоты = 1
		;mov dword ptr DS:[380h], 123B64Fh		;Timer Initial Count: 1,3 GHz / 68Hz
		mov dword ptr DS:[380h], 30000000		;new frequency for interview result showing
		smov DS, DDATA
		mov dword ptr DS:[0FCh], 0

		pop cx
		pop DS
		ret
init_LAPIC endp

spur_intr:  
iret  ;Ends without EndOfInterrupt! 

I20		equ $-_TEXT32START
i20_localtimer_handler:	
		push DS
		push ax
		smov DS, DLAPIC
		mov dword ptr DS:[0B0h], 0
		smov ds, DDATA
		inc dword ptr DS:[0FCh]
		pop ax
		pop DS
iretd


;========================================================================
;---------------------------------init_INTERVIEW------------------------
;========================================================================
rebuf proc near
		cmp di, EBUF
		jnz exit
		mov di, SBUF
exit:	ret
rebuf endp
;------------------------------------------------------------------------
pushbuf proc near
		push di
		push es
		smov es, DDATA, <_treg equ bx>	;вполнить смов используя не ax (который по умолчанию, но сейчас занят) , а bx
		mov di, word ptr es:[EBUFADR]
		stosb
		call rebuf
		cmp di, word ptr es:[SBUFADR]
		jz exit1
		mov word ptr es:[EBUFADR], di
	exit1:
		pop es
		pop di
		ret
pushbuf endp
;-----------------------------------------------------------------------
request proc near
		mov cx, 0FFF0h
w0:	
		in al, 064h
		test 2, al
		jz w1
		loop w0
		;call stop
w1:
		move_out dx, bl
		ret
request endp
;------------------------------------------------------------------------
respond proc near
		mov cx, 0FFFFh
		;mov ecx, 0FFFFFFh
w2:
		in al, 064h
		test 1, al
		jnz w3
		dec	ecx
		jnz w2
		;call stop
w3:
		in al, 060h
		ret
respond endp
;------------------------------------------------------------------------

init_INTERVIEW proc near
		push	DS
		push	ax
		push	dx		
		push	cx		
		push	bx		
		push	si		

;-----------------K E Y B O A R D - I N T E R V I E W-------------------
;keyboard_interview:

		call cout_kb
		
		req_resp 64h, 0AAh, 55h				;сброс и автотест 8042
		req_resp 64h, 0ABh, 0				;тест синхронизации с 8031
		req_resp 64h, 0AEh					;разрешить клавиатуру
		req_resp 60h, 0FFh, ACK, TESTOK		;сброс и тест 8031

;http://www.win.tue.nl/~aeb/linux/kbd/scancodes-10.html#keyboardid
;An XT keyboard does not reply
;an AT keyboard only replies with an "ACK".
;An MF2 AT keyboard reports ID "ab 83". Translation turns this into "ab 41".

		req_resp 60h, 0F2h, 0FAh			;запрос 2х байтного ID
		call respond						;got 1st byte, lets check it

		cmp eax, 000200abh					;ab (MF)
		je	its_mf
				cmp eax, 000200fah				;fa==ACK (AT)
				je	its_at
					;otherwise its XT (if it doesnt reply)
					call cout_xt
					call respond	;just in case lets get 2nd byte of id (but we dont really need it)
					jmp mouse_i

				its_at:	
					call cout_at	
					call respond	;just in case lets get 2nd byte of id (but we dont really need it)
					jmp	mouse_i

		its_mf:
			call cout_mf
			call respond		;lets get 2nd byte and check it

;------------------M O U S E - I N T E R V I E W-------------------
mouse_i:
		mov	al, ';'
		call pushbuf

		call cout_ms
		
		req_resp 64h, 0A8h						;разрешить мышь
		req_resp 64h, 60h						;запись управл. регистра
		req_resp 60h, 0							;запрещаем IRQ1, IRQ12  разрешение интерфейсов клав. и мыши
		req_resp 64h, 0D4h						;передать след. байт мыши
		req_resp 60h, 0FFh, ACK, TESTOK, 0		
		req_resp 64h, 0D4h						;передать след. байт мыши

;http://www.win.tue.nl/~aeb/linux/kbd/scancodes-13.html
;Read mouse ID command is answered by 00.
;In Intellimouse mode, the Read mouse ID command returns 03.
;In Intellimouse Explorer mode, the Read mouse ID command returns 04.

		req_resp 60h, 0F2h, 0FAh				;получить ID (ответы FA (ack) и 00 (id))
		call respond
		
		cmp	al, 0				;ps/2 id at initialization = 0
		je	its_ps2
				cmp al, 3				;imps id = 3
				je	its_imps1
						cmp al, 4				;imps2 id = 4
						je	its_imps2
						jmp uart_i
						
						its_imps2:
							call cout_im2
							jmp	uart_i

				its_imps1:
					call cout_im1
					jmp uart_i

		its_ps2:
			call cout_ps2

;--------------------U A R T - I N T E R V I E W-------------------
uart_i:
		mov al, ';'
		call pushbuf
		call cout_uart

;Test UART chip type:

		send_com <COM4+FCR>, 0E7h
		mov		dx, COM4+IIR
		in		al, dx
		test	al, 64d						;check 6th bit
		jz oldtype
				test	al, 128d				;check 7th bit
				jz t16550
						test	al, 32d				;check 5th bit
						jz t16550A
								call cout_type16750
								jmp	@f

						t16550A:
							call cout_type16550A
							jmp @f
				t16550:
					call cout_type16550
					jmp @f

	oldtype:			;Else we know the chip doesnt use FIFO, so we need to check the scratch register
		send_com <COM4+SR>, 05Ah			;Arbitary value (can be any but not 0xFF or 0x00)
		mov		dx, COM4+SR
		in		al, dx
		cmp		al, 05Ah			;If the arbitrary value comes back identical, UART is 16450
		je t16450
				call cout_type8250
				jmp @f
		t16450:	
			call cout_type16450

;Test UART speed:
@@:		
		send_com <COM4+LCR>, 128d			;set DLAB
		send_com <COM4+DLL>, 12d			;check 9600 bps
		send_com <COM4+DLH>, 0
		mov		dx, COM4+DLL
		in		al, dx
		cmp		al, 12d
		jne s2400 

				send_com <COM4+DLL>, 6d			;check 19200 bps
				send_com <COM4+DLH>, 0
				in		al, dx
				cmp		al, 6d
				jne s9600

						send_com <COM4+DLL>, 1		;check 115200 bps (top speed)
						send_com <COM4+DLH>, 0
						in		al, dx
						cmp		al, 1
						jne s19200

								call cout_115200bps
								jmp @f

						s19200:
							call cout_19200bps
							jmp @f
				s9600:
					call cout_9600bps
					jmp @f
		s2400:
			call cout_2400bps

@@:		send_com <COM4+LCR>, 0

;----------------------------------------------
		mov si, 10d
@@:		call space				;fill free space in buffer with spaces
		dec si
		jnz @b

		pop	si
		pop	bx
		pop	cx
		pop	dx
		pop	ax
		pop	DS
		ret
init_INTERVIEW endp



;========================================================================
;---------------------------------init_UART-----------------------------
;========================================================================
;При настройке UART необходимо задать: размер слова (5..8 бит); 
;									   число т.н. «стоп-битов» (1 или 1.5(2) бита), 
;									   способ задания бита чётности, 
;									   скорость передачи данных (от 50 до 115200 бод). 
ComRead proc near
		xor ax, ax
		xor bx, bx
		
		mov dx, COM1 + LSR
		mov ecx, 0FFFFFFh
w4:	
		in	al, dx
		and al, 8Fh		
		jnz	w5
		dec ecx
		jnz w4
w5:
		and al, 08Eh
		jnz stop
		mov dx, COM1 + RBR
		in al, dx

		move_out 0E9h, al
		cmp al, 'M'
		jne skppp

		mov al, '!'
		move_out 0E9h, al
skppp:
		ret
ComRead endp
;------------------------------------------------------------------------
ComWrite	proc near

		mov dx, COM1 + LSR
		mov cx, 0FFFFh
xxx:	in      al, dx
		and		al, 0AEh
		jnz		xxy
		loop xxx
xxy:
		and		al, 08Eh
		jnz		stop

		mov dx, COM1 + THR
		mov al, 'M'
		out dx, al

	    ret
ComWrite endp
;------------------------------------------------------------------------
sendcomport proc near
		push		dx
		add			bx, cx		;its implied we have COM address in 'cx' and it's REGISTER in bx(bl)
		mov			dx, bx
		out			dx, al		;in 'al' we have value
		pop			dx
		ret
sendcomport endp

;Биты 2 и 1 регистра FCR используются при первоначальной настройке порта, 
;чтобы очистить содержимое FIFO, если туда что-то уже успело попасть; 
;это вполне может случиться, так как включение устройства на другой стороне линии 
;могло выполниться гораздо раньше и оно уже могло пытаться передать какие-либо данные..

init_UART proc near
		push	si
		push	ax
		push	bx
		push	cx
		push	DS

		xor		si, si

;http://www.sci.muni.cz/docs/pc/serport.txt
unicorns:
		mov		cl,	byte ptr [COM_TABLE+si]
		inc		si
		mov		ch, byte ptr [COM_TABLE+si]
		inc		si
		push	si

		xor		si, si

		comloop:
					mov		al, [COM_INIT+si]
					inc		si
					cmp		al, 0FEh
					jz		endcom

					mov		bl, [COM_INIT+si]			
					inc		si
					call	sendcomport
					jmp		comloop
		endcom:

		pop		si
		cmp		si, 8
		jnz unicorns
	
		send_com <COM1+MCR>, 3			;set 0,1 bits of MCR (DTR, RTS)
		send_com <COM2+MCR>, 3			;set 0,1 bits of MCR (DTR, RTS)
		send_com <COM3+MCR>, 3			;set 0,1 bits of MCR (DTR, RTS)
		send_com <COM4+MCR>, 3			;set 0,1 bits of MCR (DTR, RTS)

		;call ComRead

		pop		DS
		pop		cx
		pop		bx
		pop		ax
		pop		si
		ret
init_UART endp


;========================================================================
;---------------------------------FINAL_CYCLE--------------------------- by IVAN, (c)
;========================================================================
FINAL_CYCLE proc near
		push ax
		push edi
		push es

		smov		ds, DDATA
		mov			word ptr DS:[SBUFADR], SBUF		;start buffer
		mov			word ptr DS:[EBUFADR], SBUF
		
inf_cycle:
		mov			dword ptr ds:[0FCh], 0
		;		xchg bx, bx
		call		init_INTERVIEW
		;		xchg bx, bx
aw0:
		
		cmp			dword ptr ds:[0FCh], 0		;connect interviews with timer's freq
		jz		aw0

		mov			di, ds:[SBUFADR]	
aw1:
		cmp			di, word ptr ds:[EBUFADR]		;check if need to restart buffer
		jnz		aw2

		move_out	0E9h, 10						;\n
		move_out	0E9h, 13
		jz		inf_cycle
aw2:
		mov			al, ds:[di]
		inc			di
		call		rebuf
		mov			word ptr ds:[SBUFADR], di
		move_out	0E9h, al
		jmp		aw1									;come back with renewed buffer

		ret
FINAL_CYCLE endp


;========================================================================
;---------------------------------stop----------------------------------
;========================================================================
stop    proc near
        cli
        hlt
        jmp     short stop
stop    endp

end32	equ	$-_TEXT32START
_TEXT32 ends



;----------------------------------------------------------------------
;чтобы было понятно что все начинается именно здесь а !не! на верху в старте32

;__$$$$____$$__$$__$$$$$__$$$$$__$$_____$$_$$____$$_ 
;_$$_______$$__$$_$$___$$_$$___$_$$_____$$_$$$__$$$_ 
;_$$_______$$$$$$_$$___$$_$$$$$__$$$$$__$$_$$_$$_$$_ 
;_$$_______$$__$$_$$___$$_$$___$_$$___$_$$_$$____$$_ 
;__$$$$____$$__$$__$$$$$__$$$$$__$$$$$__$$_$$____$$_ 
;___________________________________________________ 
;__$$$$$$___$$$$$_____$$$$$$_____$$$$$____$$____$$__ 
;__$$______$$___$$____$$__$$____$$___$$___$$$__$$$__ 
;__$$______$$___$$____$$__$$____$$___$$___$$_$$_$$__ 
;__$$______$$___$$____$$__$$____$$___$$___$$____$$__ 
;__$$_______$$$$$_____$$$$$$_____$$$$$____$$____$$__ 
;___________________$$______$$______________________ 



_TEXT   segment byte public 'CODE' use16
assume cs:_TEXT, ds:nothing
_TEXTSTART = $
start:
		cli
		lss		SP, dword ptr STKPTR
		sti
		call init_BIOS

		jmp protected_mode_on


;========================================================================
;------------------------------init_BIOS--------------------------------
;========================================================================


STKPTR      dw		0FFFEh, 09000h
init_BIOS:  mov  	dx, 0C000h	
cycle: 		mov  	DS, dx
  			mov  	ax, 80h			
  			cmp  	word ptr DS:[0], 0AA55h	
  			jnz  	nxt	
			call scanbios

  			movzx  	ax, byte ptr DS:[2]
  			add  	al, 3h		
  			and  	al, 0FCh
			shl  	ax, 5  	

nxt: 		add  	dx, ax			
  			cmp  	dx, 0F000h			
  			jb  	cycle		
		ret
  
scanbios proc near
  		cld	
  		xor     si, si
 		xor     cx, cx
  		mov     ch, DS:[2]		
  		xor     bl, bl
chcksm: lodsw	
        add     al, ah			
        add     bl, al			
        dec     cx
        jnz     short chcksm
        or      bl, bl	
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


;========================================================================
;-----------------------------DESCRIPTORS-------------------------------
;========================================================================
;
; Usage:
;	.descriptor								# all zeroes descriptor, must be #0 in table
;	.descriptor	limit=0x12345678, base=0x9ABCDEF0, r=1, p=1, dpl=0	# flag 'r' used, 'w' is'nt used => code segment
;	.descriptor	limit=0x12345678, base=0x9ABCDEF0, w=1, p=1, dpl=0	# flag 'r' is'nt used, 'w' is used => data segment
;
; g - granularity (0 or 1), set to 1 automatically, when limit > 0x100000; elsewhere g defines bytes/pages for limit field
; x - default operation size (0 or 1); means 32-bit instructions for code segment 
														;or 32-bit pointers for stack segment 
														;or 2^16 vs 2^32 limit for expand down segments
; l - 64-bit long mode segment
; p - presence in memory (0 or 1)
; dpl - DPL
; a - accessed (0 or 1)
; only for code segments:
;  r - execute only or execute/read (0 or 1) NOTE: r=0 or r=1 assumes code segment
;  c - (only for code segment, 0 or 1) conforming segment
; only for data segments:
;  w - read only or read/write (0 or 1) NOTE: w=0 or w=1 assumes data segment
;  ed - (only for data segment, 0 or 1) expand down segment
; type - 1,3,9,11 or 2 can be used only (TSSes, LDT)
;  type, r and w bits are mutually defined
;бит гранулярности меняет лимит и он уже не 0Fh
GDT:
		dd 0, 0
		descriptor _limit=0Fh, _base=0F0000h, _g=1, _x=1, _r_w=1, _d_c=1 ;32х разрядный сегмент кода
		descriptor _limit=100h,	_base=400h, _r_w=1, _g=0				;16х данные BIOS (BDA)
		descriptor _limit=0FFFFh, _base=90000h, _g=0					;32х разрядный стек
		descriptor _base=0FEC00000h, _g=1, _x=1
		descriptor _base=0FEE00000h, _g=1, _x=1
GDT_SIZE equ $-GDT
GDTR:	dw GDT_SIZE-1
		dd STARTGDT SHL 4

;----------------------------------------------------------------------
IDT:
		intdescriptor _offset=I20, _selector=8, _use32=1
IDT_SIZE equ $-IDT+256
IDTR:	dw IDT_SIZE-1
		dd STARTIDT SHL 4


;========================================================================
;--------------------------PROTECTED-MODE-ON----------------------------
;========================================================================
;Изменить первую работу этого семестра (N1) так, чтобы:
;a) Процессор был переведён в режим защищённого сегментного преобразования.
;b) Финальный цикл выполнялся в 32-x битовом со смещениями не более 216 режиме.
;c) Обработка прерываний осуществлялась в 32-x битовом со смещениями не более 216 режиме.
;d) Стек размещался в 32-x битовом со смещениями не более 216 обычном сегменте.
;e) Данные биос (BDA, Bios Data Area) размещались в 16-ти битовом сегменте.
protected_mode_on:
		call open_A20			;Для обращенний за пределы первого мегабайта надо включить режим «Линия А20 разрешена»
		call copy_gdt
		call copy_idt
		call restrict_NMI
		
		mov	eax, cr0			;Очищаем бит PE в CR0 (выключаем сегментацию)
		or	eax, 1
		mov	cr0, eax

		db 0EAh
		dw offset start32, 8		;прыжок в 32битный код наверх на старт32

;----------------------------------------------------------------------
copy_gdt:
		smov ds, cs
		mov si, offset GDT
		smov es, STARTGDT
		xor di, di
		mov cx, GDT_SIZE/4
		cld
		rep movsd
		lgdt fword ptr GDTR
		ret

;----------------------------------------------------------------------
copy_idt:
		mov si, offset IDT
		smov es, STARTIDT
		mov di, 256			;пропускаем 21 интеррапт и копируем только интеррапт таймера лапика
		mov cx, IDT_SIZE/4	;cx=2
		cld
		rep movsd
		lidt fword ptr IDTR
		ret

;----------------------------------------------------------------------
restrict_NMI:
		cli
		in al, 70h
		or al, 80h
		out 70h, al
		ret

;----------------------------------------------------------------------
open_A20:
		in al, 92h
		or al, 2
		out 92h, al
		ret



;------------------------------------------------------------------------
org	0FFF0h-end32
	db	0EAh					; JMP FAR
	dw	offset start			; offset
	dw	0F000h+(end32 SHR 4)	; segment

org	0FFFEh-end32
	dw	99FCh					; PC 

_TEXT	ends
end		start

