.model flat
.code

;Add256BE (c, a, b)
;c = a + b
;typedef int256_x int[8]
;Big-Endian

public func
func proc c public
	push ebp
	mov 		ebp, esp
	push esi
	push ebx
			
	mov edx, 8[ebp]		;c
	mov ebx, 12[ebp]	;a
	mov esi, 16[ebp]	;b

;	mov ecx, 6
;loop:	mov eax, [ebx+ecx*04h]
;	adc eax, [esi+ecx*04h]
;	mov [edx+ecx*04h], eax
;	каждый раз декремент/пуш	

	mov eax, 28[ebx]	;eax=a[7]
	add eax, 28[esi]	;+b[7]
	mov 28[edx], eax	;c[7]=a[7]+b[7]
	
	mov eax, 24[ebx]	;eax=a[6]
	adc eax, 24[esi]	;+b[6]
	mov 24[edx], eax	;c[6]=a[6]+b[6]

	mov eax, 20[ebx]	;5
	adc eax, 20[esi]	
	mov 20[edx], eax

	mov eax, 16[ebx]	;4
	adc eax, 16[esi]	
	mov 16[edx], eax

	mov eax, 12[ebx]	;3
	adc eax, 12[esi]	
	mov 12[edx], eax

	mov eax, 8[ebx]		;2
	adc eax, 8[esi]	
	mov 8[edx], eax

	mov eax, 4[ebx]		;1
	adc eax, 4[esi]	
	mov 4[edx], eax

	mov eax, [ebx]		;0
	adc eax, [esi]	
	mov [edx], eax	

	pop ebx
	pop esi
	mov eax, 16[esp]
	pop ebp
	ret
func endp
end

