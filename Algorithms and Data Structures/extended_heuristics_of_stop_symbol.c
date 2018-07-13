#include <stdio.h>

int max(int a,int b) {
    if (a>b) {
        return a;
    }
    else
        return b;
}

void delta1(char* s,int *d) {
}

void delta(char* s,int** d) {
}

int main(int arg,char** argv)
{
    char s[1000],t[1000];
    if (arg==3) {
        strncpy(s,argv[1],999);
        strncpy(t,argv[2],999);
    }
    else {
        gets(s);
        gets(t);
    }
    int n=strlen(s);
    int m=strlen(t);
    int k=n-1,i;
    int d[n][256];
    /* 
        for (j=0;j<=k;j++) {
            d[i][s[i+j]]=j;
        }
    }*/
    
    for (i=0;i<n;i++) {
        int j;
        int k=n-i;
        for (j=0;j<256;j++) {
            d[i][j]=k;
        }    
    }
    int j;
    for (i=0;i<256;i++){
        int k=n;
        int j;
        for (j=0;j<n;j++){
            if (i==s[j]) {
                if (k!=n) 
                    d[j][i]=n-k-1;
                k=j;
                continue;
            }
            if (k!=n) 
                d[j][i]=n-k-1;
        }
    }
    
    while (k<m) {
        i=n-1;
        while (t[k]==s[i]) {
            if (i==0) {
                printf("%d",k);
                return 0;
            }
            i--;
            k--;
        }

        k+=d[i][t[k]];
    }
    printf("%d",m);
    return 0;
}