#include <stdio.h>
 
int main(void)
{
  int n,i,j,max=0,k;
  scanf("%d",&n);
  n=abs(n);
  k=(int)(sqrt(n))+2;
  int A[k];
  for (i=2;i<=k;i++) {
    A[i]=1;
  }
  for (i=2;i<=(int)(sqrt(k))+1;i++) {
    if (A[i]==1)
      for (j=i*i;j<=k;j+=i) {
        A[j]=0;
      }
  }
  for (i=2;i<=k;i++) {
    if (A[i]==1) {
      while (!(n%i)) {
        max=i;
        n=n/i;
      }
    }
  }
  if (n>1) 
    max=n;
  printf("%d",max);
  return 0;
}




#include <stdio.h>
#include <math.h>

int main()
{
        int x;
        scanf ("%d", &x);
        int t=ceil(sqrt(abs(x))+1);
        int array[t];
        int i, j;
        for (i=0; i<=t; i++) {
                array[i]=i;
        }
        for (i=2; i<=t; i=i++) {
                if(array[i]!=0) {
                        for (j=i*2; j<=t; j=j+i) {
                                array[j]=0;
                        }
                }
        }
        array[1]=0;
        for (i=2; i<t; i++) {
                if (((array[i])!=0) && (x%array[i]==0) && (x!=array[i]))
                        x=x/array[i];
                if ((x % i == 0) && (x != array[i])) {
                        while (x % i == 0)
                        x = x/array[i];
                }
        }
        x=abs(x);
        printf ("%d", x);
        return 0;
}