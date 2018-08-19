.model flat
.code

;EXAM(17june)
;MAX256LE (c, a, b)
;c = max( a, b )
;typedef int int256_x[8]
;LE

public func
func	proc c public
	push		ebp
	mov		ebp, esp
	push		esi
	push		ebx
	push		ecx
	
	mov		esi, 8[ebp]
	mov		edx, 12[ebp]
	mov		ebx, 16[ebp]

	mov		eax, [edx]
	cmp		eax, [ebx]
	jl	loopb
	jne	loopa

	mov		ecx, 7
comp:
	mov		eax, [edx+ecx*4]
	cmp		eax, [ebx+ecx*4]
	ja	loopa
	jb	loopb
	dec		ecx
	jnl	comp

loopa:	mov		ecx, 7
a:	mov		eax, [edx+ecx*4]
	mov		[esi+ecx*04h], eax
	dec		ecx
	jnl	a
	jmp	exit

loopb:	mov		ecx, 7			
b:	mov		eax, [ebx+ecx*4]
	mov		[esi+ecx*04h], eax
	dec		ecx
	jnl	b

exit:
	pop		ecx
	pop 		ebx
	pop 		esi
	pop		ebp
	
	ret

func	endp
end


//Main.c:

#include <stdio.h>
#pragma warning(disable:4100)

typedef int int256[8];

int main( int ac, char **av, char **env )
{
	int256   a, b, c;
	a[0] = 0x00000000; a[1] = 0x00000000; a[2] = 0x00000000; a[3] = 0x00000000; a[4] = 0x00000000; a[5] = 0x00000000; a[6] = 0x00000000; a[7] = 0x00000000;
	b[0] = 0x00000000; b[1] = 0xF0000000; b[2] = 0x00000000; b[3] = 0x00000000; b[4] = 0x00000000; b[5] = 0x00000000; b[6] = 0x00000000; b[7] = 0x00000000;
	c[0] = 0x77777777; c[1] = 0x77777777; c[2] = 0x77777777; c[3] = 0x77777777; c[4] = 0x77777777; c[5] = 0x77777777; c[6] = 0x77777777; c[7] = 0x77777777;
	func(c, a, b);
	printf( "0x%08X%08X%08X%08X%08X%08X%08X%08X, \n", a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7] );
	printf( "0x%08X%08X%08X%08X%08X%08X%08X%08X \n", b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7] );
	printf( "0x%08X%08X%08X%08X%08X%08X%08X%08X\r\n\n", c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7] );


	a[0] = 0xF0000000; a[1] = 0x00000000; a[2] = 0x00000000; a[3] = 0x00000000; a[4] = 0x00000000; a[5] = 0x00000000; a[6] = 0x00000000; a[7] = 0x00000000;
	b[0] = 0x00000000; b[1] = 0xF0000000; b[2] = 0x00000000; b[3] = 0x00000000; b[4] = 0x00000000; b[5] = 0x00000000; b[6] = 0x00000000; b[7] = 0x00000000;
	c[0] = 0x77777777; c[1] = 0x77777777; c[2] = 0x77777777; c[3] = 0x77777777; c[4] = 0x77777777; c[5] = 0x77777777; c[6] = 0x77777777; c[7] = 0x77777777;
	func(c, a, b);
	printf( "0x%08X%08X%08X%08X%08X%08X%08X%08X, \n", a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7] );
	printf( "0x%08X%08X%08X%08X%08X%08X%08X%08X \n", b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7] );
	printf( "0x%08X%08X%08X%08X%08X%08X%08X%08X\r\n\n", c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7] );


	a[0] = 0x00000000; a[1] = 0x00000000; a[2] = 0x00000000; a[3] = 0x00000000; a[4] = 0x00000000; a[5] = 0x00000000; a[6] = 0x00000000; a[7] = 0x00000000;
	b[0] = 0xF000000; b[1] = 0x00000000; b[2] = 0x00000000; b[3] = 0x00000000; b[4] = 0x00000000; b[5] = 0x00000000; b[6] = 0x00000000; b[7] = 0x00000000;
	c[0] = 0x77777777; c[1] = 0x77777777; c[2] = 0x77777777; c[3] = 0x77777777; c[4] = 0x77777777; c[5] = 0x77777777; c[6] = 0x77777777; c[7] = 0x77777777;
	func(c, a, b);
	printf( "0x%08X%08X%08X%08X%08X%08X%08X%08X, \n", a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7] );
	printf( "0x%08X%08X%08X%08X%08X%08X%08X%08X \n", b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7] );
	printf( "0x%08X%08X%08X%08X%08X%08X%08X%08X\r\n\n", c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7] );


	a[0] = 0x00000000; a[1] = 0x00000000; a[2] = 0x00000000; a[3] = 0x00000000; a[4] = 0x00000000; a[5] = 0x00000000; a[6] = 0x00000000; a[7] = 0x00000000;
	b[0] = 0xF0000000; b[1] = 0x00000000; b[2] = 0x00000000; b[3] = 0x00000000; b[4] = 0x00000000; b[5] = 0x00000000; b[6] = 0x00000000; b[7] = 0x00000000;
	c[0] = 0x77777777; c[1] = 0x77777777; c[2] = 0x77777777; c[3] = 0x77777777; c[4] = 0x77777777; c[5] = 0x77777777; c[6] = 0x77777777; c[7] = 0x77777777;
	func(c, a, b);
	printf( "0x%08X%08X%08X%08X%08X%08X%08X%08X, \n", a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7] );
	printf( "0x%08X%08X%08X%08X%08X%08X%08X%08X \n", b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7] );
	printf( "0x%08X%08X%08X%08X%08X%08X%08X%08X\r\n\n", c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7] );

	return 0;
}

