#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
char* fibstr(int n)
{
        char* first_str;
        char* prev;
        char* current;
        int i = 0, next_string_len;
 
        prev = malloc(1000);
        strcpy (prev,"a");
 
        current = malloc(1000);
        strcpy (current,"b");
 
        if (n == 1)
        return prev;
        else
        if (n == 2)
            return current;
        else {
             for (i=3;i<=n;i++){
             char* next_str;
 
             next_string_len = strlen(prev) + strlen(current) + 1;//
             next_str = (char* )malloc(next_string_len);  
                     memset(next_str,0,next_string_len);     
 
                     strcat(next_str,prev);       
                     strcat(next_str,current);   
 
                     free (prev);         
 
             prev = current; 
             current = next_str;
             }
        }
        free (prev); 
 
        return current;
}
 
int main()
{
        unsigned long n;
        scanf("%d", &n);
 
        char* result;
        result = fibstr(n);
        printf("%s", result);
 
        free(result);
        return 0;
}