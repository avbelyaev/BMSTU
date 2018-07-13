#include <stdio.h>

int main(void)
{
  int n,m;
  scanf("%d %d",&n,&m);
  int A[n][m];
  int Ymin[m],Xmax[n];
  int i,j;
  for (i=0;i<=m-1;i++) {
    Ymin[i]=0;
  }
  for (i=0;i<=n-1;i++) {
    Xmax[i]=0;
  }
  for (i=0;i<=n-1;i++) {
    for (j=0;j<=m-1;j++) {
      scanf("%d",&A[i][j]);
      if (A[i][j]>A[i][Xmax[i]])
        Xmax[i]=j;
      if (A[i][j]<A[Ymin[j]][j])
        Ymin[j]=i;
    }
  }
  for (i=0;i<=n;i++) {
    if (i==n) 
      printf("none");
    else
      if (Xmax[Ymin[i]]==i) {
        printf("%d %d",Ymin[i],i);
        break;
    }
  }
  return 0;
}