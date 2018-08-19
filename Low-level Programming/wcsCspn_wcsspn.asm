#include <stdio.h>
#include <stdlib.h>
#include <wchar.h>
#include <string.h>


size_t x_wcscspn(char *s, char *reject)
{
	if (!s || !(*s) || !reject || !(*reject)) return strlen(s);
	int len = strlen(reject), len2 = strlen(s);
	char* temp;
	char* one = s;
	int templen;
	do {
		__asm__ __volatile__ ("repne scasb \n\t" : "=c"(templen), "=D"(temp) : "a"(*s), "D"(reject), "c"(len));
		__asm__ __volatile__  goto("jz %l[s_ret] \n\t" : : : "cc" : s_ret);
	} while (*(++s));
    s = 0;
s_ret:
    if (s == 0) return len2;
	return s - one;
}

size_t x_wcsspn(wchar_t *s, wchar_t *accept)
{
	if (!s || !(*s) || !accept || !(*accept)) return wcslen(s);
	int templen, len = wcslen(accept), len2 = wcslen(s);
	wchar_t* temp;
	wchar_t* start = s;
	//this works on Win only. in Linux version replace "scasw" with "scasl"
	do {    __asm__ __volatile__ ("repne scasw \n\t" : "=c"(templen), "=D"(temp) : "a"(*s), "D"(accept), "c"(len));
            __asm__ __volatile__  goto("jnz %l[exit] \n\t" : : : "cc" : exit);
	} while (*(++s));

    s = 0;
    exit:
    if (s == 0) return len2;
	return (s - start);
}

/*asm [volatile] ( AssemblerTemplate
                      : OutputOperands
                      [ : InputOperands
                      [ : Clobbers ] ])

     asm [volatile] goto ( AssemblerTemplate
                           :
                           : InputOperands
                           : Clobbers
                           : GotoLabels)*/

//"cc"
//  The "cc" clobber indicates that the assembler code modifies the flags register. On some machines, GCC represents the condition codes as a specific hardware register; "cc" serves to name this register. On other machines, condition code handling is different, and specifying "cc" has no effect. But it is valid no matter what the target.

int main()
{
    char str[10] = "abcddddef";
    wchar_t wstr[10] = L"abcddddef";
    char rej[10] = "efhg";
    wchar_t acc[10] = L"cabd";
    int ans;
    ans = x_wcscspn(str, rej);      //CHAR
    printf("wcsCspn=%d\n", ans);
    ans = x_wcsspn(wstr, acc);      //WCHAR_T
    printf("wcsspn=%d\n", ans);
    return 0;
}

