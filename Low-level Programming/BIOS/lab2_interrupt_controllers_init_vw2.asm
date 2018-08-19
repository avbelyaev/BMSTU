;Инициализировать контроллеры прерываний в Virtual Wire II (с IOAPIC) режиме.
;Настроить LocalTimer для генерации прерываний с вектором 0x20 на частоте 68Гц.

.586
_TEXT   segment byte public 'CODE' use16
assume cs:_TEXT, ds:nothing
org     100h			
start:
        cli		
        lss     SP, dword ptr STKPTR
        sti		

	call  	init_BIOS
	call	init_PIC
	call	unreal_mode_on
	call	init_IOAPIC
	call	init_LAPIC
	call 	stop

;========================================================================
;---------------------------------init_BIOS-----------------------------
;========================================================================

STKPTR  dw	0FFFEh, 09000h

init_BIOS:
	mov  	dx, 0C000h
		
cycle: 	mov  	DS, dx
  	mov  	ax, 80h			
  	cmp  	word ptr DS:[0], 0AA55h	
  	jnz  	nxt	
		
	call scanbios

  	movzx  	ax, byte ptr DS:[2]
  	add  	al, 3h		
  	and  	al, 0FCh
	shl  	ax, 5  	
		 
nxt: 	add  	dx, ax			
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
;---------------------------------init_PIC-------------------------------
;========================================================================

init_PIC proc near

;1) начальный сброс (безличный EOI)(запись OCW2 с кодом безличного EOI — 0x20 в чётные порты обоих контроллеров)

	mov al, 20h		;pic0
	out 20h, al
	mov al, 20h		;pic1
	out 0A0h, al

;2) Загрузка ICW1..ICW4 в оба контроллера
; ICW1 записывается в чётный порт контроллера, после чего в нечетный порт должны быть немедленно записаны ICW2..ICW4. 
; ICW3 надо указывать только если используется каскадирование (так и есть), а ICW4 только если бит ICW4 в ICW1 установлен.

	mov al, 11h		;pic0
	out 20h, al
	mov al, 11h		;pic1
	out 0A0h, al

; ICW2 записывается в нечётный порт контроллера и задает номер вектора прерывания ЦПУ, соответствующий нулевой линии запроса прерывания. 
; Линии (1..7) будут отображены на вектора N+1, .., N+7. В реальном режиме это обычно вектора 0x08 для PIC #0 и 0x70 для PIC #1.

	mov al, 08h		;pic0
	out 21h, al
	mov al, 70h		;pic1
	out 0A1h, al
; ICW3 записывается в нечётный порт контроллера и задает:
; для PIC #0 — битовую маску линии запроса, к которой подключен PIC #1 (IRQ2, маска 0b00000100)
; для PIC #1 — номер уровня ведомого контроллера (обычно 2)

	mov al, 04h		;pic0
	out 21h, al
	mov al, 02h		;pic1
	out 0A1h, al

;ICW4 записывается в нечётный порт контроллера, если бит ICW4 в ICW1 был установлен.

	mov al, 05h		;pic0
	out 21h, al
	mov al, 01h		;pic1
	out 0A1h, al

;3) разрешение выбранных IRQ
; OCW1 чтение/запись в нечётный порт контроллера и задает битовую маску запрещенных линий IRQ. Обычно запрещены линии IRQ 3,4,5,7 на PIC #0 (маска 0xB8) и 8,9,A,B,F на PIC #1 (маска 0x8F)
;PIC0
	mov al, 0FBh 		;1111 1011
	out 21h, al
;PIC1
	mov al, 0FFh		;запретили все
	out 0A1h, al
ret
init_PIC endp

;========================================================================
;---------------------------------unreal_mode_on-------------------------
;========================================================================
.586p
unreal_mode_on proc near
		push	ds
		push	es
		push	gs
		push	fs
		push	cx
		push	ax
		push	si
		push	di

		mov	ax, cs
		mov	ds, ax
		mov	si, offset _TEXT:unreal_mode_c

		mov	ax, 9000h
		mov	es, ax
		xor	di, di

		cld
		mov	cx, 6
		rep	movsd
		lgdt fword ptr _TEXT:unreal_mode_b

		mov	eax, cr0
		or	eax, 1
		mov	cr0, eax
		db	0EAh
		dw	$+4, 8h

		mov	ax, 10h
		mov	ss, ax
		mov	ds, ax
		mov	es, ax
		mov	fs, ax
		mov	gs, ax

		mov	eax, cr0
		and	eax, 0FFFFFFFEh
		mov	cr0, eax
		db	0EAh
		dw	$+4, 0F000h

		mov	ax, 9000h
		mov	ss, ax

		in	al, 92h
		test	al, 2
		jne	unreal_mode_d
		or	al, 2
		jmp short $+2
		out	092h, al
unreal_mode_d:
		lgdt fword ptr _TEXT:unreal_mode_a

		pop	di
		pop	si
		pop	ax
		pop	cx
		pop	fs
		pop	gs
		pop	es
		pop	ds
		ret	0
unreal_mode_on endp

unreal_mode_a:
	dw	0
	dd	0, 0

unreal_mode_b:
	dw	23
	dd	90000h, 0

unreal_mode_c:
	dd	0, 0, 0FFFFh, 8F9A0Fh, 0FFFFh, 8F9200h

;========================================================================
;---------------------------------init_IOAPIC----------------------------
;========================================================================

;1.IOAPIC обычно отображается на адрес IOAPICBASE 0xFEC00000 (адрес может быть несколько изменен средствами северного моста)
;2.IOAPIC представлен двумя регистрами: "селектор" (selector) по адресу IOAPICBASE и "окно" (window) по адресу IOAPICBASE+0x10
;3.При работе с IOAPIC в селектор записывается номер регистра IOAPIC, а из окна производится чтение или запись регистра
;4.Термин "Reserved", определяющий те или иные биты регистров IOAPIC обозначает не "MBZ", а то, что содержимое этих разрядов нельзя изменять, 
;при необходимости изменить регистр надо сначала прочитать прежнее значение, потом только нужные биты обнулить или установить, а затем записать обратно.
.586p
init_IOAPIC proc near

	push	DS
	push	cx

	xor	ax, ax
	mov	DS, ax

;адрес селектора:	0xFEC00000h
;адрес окна: 	 	0xFEC00010h

;24 пары 32-х разрядных регистров образуют 24 записи, по 64 бита на описание каждого IRQ; порядок little-endian (от 0x10,0x11 до 0x3E,0x3F):
;линия IOAPIC  	|   0   |   1   |  ...  |   22  |   23  |
;регистр 	| 10,11	| 12,13	|  ...	| 3C,3D	| 3E,3F	|

;линия 0 настраивается отдельно от всех, так как по ней подключен PIC. остальные линии маскируются
;0 линия, 10 рег
	mov	dword ptr DS:[0FEC00000h], 10h		;номер регистра отвечающего за 0-ю линию -> в адрес селектора (запрос) (0-я линия - 0x10,0x11 рег, 23-я линия - 0x3E,0x3F рег)
	mov 	eax, dword ptr DS:[0FEC00010h]		;работа с регистром в окне (ответ)
  	and 	eax, 0FFFE60FFh				;Физическая адресация LAPIC, обнуляем Trigger Mode(15й бит), Pin Polarity(131 бит), 
							;возбуждение прерывания передним фронтом (состояние 1) т.к. настраиваем линию irq0 - линию по которой подключен PIC
  	or 	eax, 0700h				;Delivery Mod = 111 = ExtInt - Causes the processor to respond the intr as if the intr originated in an externally connected intr controller
							;прерывания идут в неизменном виде от устройства. Если Delivery Mod = Fixed то они уже были бы заданы и обрабатывались бы заданные, а не полученные прерывания							
  	mov 	dword ptr DS:[0FEC00010h], eax		;настроили маску eax и записали ее обратно

;0 линия, 11 рег
  	mov 	dword ptr DS:[0FEC00000h], 011h		;второй регистр 0-й линии
  	mov 	eax, dword ptr DS:[0FEC00010h]	
  	and 	eax, 0FFFFFFh				;Обнулить целевой ID лапика приёмника (0-й lapic единственного процессора) т.к бочс эмулирует лишь 1 процессор
  	mov 	dword ptr DS:[0FEC00010h], eax	

;на этом 0-я линия настроена, остальные 23 линии маскируем в цикле. начинаем с 1-й линии, 12 регистра
  	mov 	ecx, 012h						
loop_masking:
  	mov 	DS:[0FEC00000h], ecx			;номер очередного регистра в селектор
  	mov 	eax, dword ptr DS:[0FEC00010h]		;работа в окне
  	or 	eax, 010000h				;маскирование (16-й бит MASK = 1)
  	mov 	dword ptr DS:[0FEC00010h], eax		;записываем обратно
  	add 	ecx, 2					;шаг цикла - через один	т.к бит маскирования(16й) находится в первом (левом) регистре линии (учесть что little-endian)
  	cmp 	ecx, 040h				;последний обработанный регистр - 3E. 3Eh + 2h = 40h			
  	jnz 	loop_masking

  	pop 	cx								
  	pop 	DS							
ret
init_IOAPIC endp

;========================================================================
;---------------------------------init_LAPIC-----------------------------
;========================================================================

init_LAPIC proc near

		push DS
		push cx

;2) Убедиться, что установлен бит разрешения LAPIC(бит номер 0x0B, маска 0x800) в регистре базового адреса LAPIC; регистр базового адреса LAPIC известен как MSR 0x1B (IA32_APPIC_BASE)

		mov	ecx, 01Bh
		rdmsr
		or	eax, 0800h
		wrmsr

		xor 	ax, ax
		mov	DS, ax

;если во время обработки какого-то прерывания появляется прерывание с повышенным приоритетом, то будет отправлен spurios vector

		mov DS:[040h], offset spur_intr			;обработка фиктивных (spurious) прерываний
		mov DS:[042h], cs				;аналогично обычным прерываниям далее

;3) Установить бит разрешения LAPIC в дескрипторе вектора фиктивных (spuriouse) прерываний

		mov dword ptr DS:[0FEE000F0h], 0110h		;8-й бит ApicEn = ApicEnabled - включение апика. 4-й бит отвечает за spurious vector (вектор фиктивных прерываний)

;4) Настроить вектора прерываний LINT0, LINT1, локального таймера и т. д.

		or dword ptr DS:[0FEE00350h], 010000h		;маскирование LINT0 
		or dword ptr DS:[0FEE00360h], 010000h		;маскирование LINT1

;обработчик обычных прерываний
		mov DS:[080h], offset interrupt_VANO		;адрес обработчика в таблице векторов прерываний: 0xVector * 0x4. смещение в сегменте cs в данном случае будет: 0x20 * 4 = 080h
		mov DS:[082h], cs				;сегмент 082h

;Local Timer: Vector 0x20, Frequency 68 Hz
		mov	eax, dword ptr DS:[0FEE00320h]		;LVT[0] - Timer
		and	eax, 06F00h				
;настройка вектора прерываний
		or	eax, 020020h				;17-й бит Timer = 1 => таймер периодический. 20h в конце маски - вектор 0x20 из задания
		mov dword ptr DS:[0FEE00320h], eax

		or dword ptr DS:[0FEE003E0h], 0Bh		;Timer Divide Configuration: 1011, d0(0й бит)=1, d1(1й бит)=1, d2(3й бит)=1 => делитель частоты = 1. 
								;т.к. частота умещается в отведенные рамки то дополнительно ее не делим (делитель = 1)
;настройка частоты таймера
		mov dword ptr DS:[0FEE00380h], 123B64Fh		;Timer Initial Count: 1,3 GHz (частота IPS из bochs) / 68 Hz (частота из задания) = 123B64Fh

		mov dword ptr DS:[012300h], 0			;обнуляем ячейку 12300h (обязательно не из 1-го мб памяти) для последующего накопления там прерываний

		pop cx
		pop DS
ret
init_LAPIC endp

spur_intr:  	;обработка вектора фиктивных прерываний
iret		;Ends without EndOfInterrupt! as said in documentation 3A, chapter 10

interrupt_VANO:	;обработка прерываний
		push DS
		push ax

		xor ax, ax
		mov DS, ax

		mov dword ptr DS:[0FEE000B0h], 0	;EOI
		inc dword ptr DS:[012300h]		;увеличиваем значение в ячейке куда "накапливаем" прерывания

		pop ax
		pop DS
iret

;========================================================================
;---------------------------------stop-----------------------------------
;========================================================================

stop	proc	near
	; cli
	hlt
	jmp	short stop
stop	endp

;--------------------------------------------------------------------------

org	0FFF0h
	db	0EAh		; JMP FAR
	dw	offset start	; offset
	dw	0F000h		; segment


org	0FFFEh
	dw	99FCh		; PC 

_TEXT	ends
end	start


;P.S. команды бочса (bochsBDG):
;"c" = continue
;Ctrl+C = stop
;Ctrl+C = fullstop
;"x 0x12300" - проверить ячейку 12300h

;P.P.S. для частоты 68Гц в ячейке 12300h за 5сек накопится значение в ~350d прерываний (15Eh)

