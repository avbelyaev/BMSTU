#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int *tree_data = 0;

#define min(x,y) (x<y?x:y)
#define max(x,y) (x>y?x:y)

#define TREE_FUNCT(x,y) (max(x,y))

void build (int a[], int v, int tl, int tr) {
        if (tl == tr)
		tree_data[v] = a[tl];
	else {
		int tm = (tl + tr) / 2;
		build (a, v*2, tl, tm);
		build (a, v*2+1, tm+1, tr);
		tree_data[v] = TREE_FUNCT(tree_data[v*2],tree_data[v*2+1]);
	}
}


int compute (int v, int tl, int tr, int l, int r, int* value) {
	int tm;
	int lvalue,rvalue;
	int lresult,rresult;

	if (l > r)
		return 0;
	if (l == tl && r == tr){
		*value = tree_data[v];
		return 1;
	}
		
	tm = (tl + tr) / 2;

	lresult = compute (v*2, tl, tm, l, min(r,tm),&lvalue);
	rresult = compute (v*2+1, tm+1, tr, max(l,tm+1), r,&rvalue);

	if ((lresult>0)&&(rresult>0)){
		*value = TREE_FUNCT(lvalue,rvalue);
	} else{
		if (lresult)
			*value = lvalue;
		if (rresult)
			*value = rvalue;
	}

	

	return 1;
}

void update (int v, int tl, int tr, int pos, int new_val) {
	if (tl == tr)
		tree_data[v] = new_val;
	else {
		int tm = (tl + tr) / 2;
		if (pos <= tm)
			update (v*2, tl, tm, pos, new_val);
		else
			update (v*2+1, tm+1, tr, pos, new_val);
		tree_data[v] = TREE_FUNCT(tree_data[v*2],tree_data[v*2+1]);
	}
}

int main(int argc,char** argv){
	int* data = 0;
	int i,N;
	int ops_count;

	scanf("%d",&N);

	data = (int*)malloc(sizeof(int) * N);
	tree_data = (int*)malloc(N * 4 * sizeof(int));


	for (i=0;i!=N;i++)
		scanf("%d",&(data[i]));

	build(data,1,0,N-1);



	scanf("%d",&ops_count);
	while(ops_count){
		char command[255];
		int l1,l2;

		memset(command,0,255);
		scanf("%s",command);
		scanf("%d",&l1);
		scanf("%d",&l2);
		
		if (command[0]=='M'){
			
			if (l2 <= N-1){
				int result;
				compute(1,0,N-1,l1,l2,&result);
				printf("%d\n",result);
			}
		}
		if (command[0]=='U'){
			if (l1 <= N-1){
				update(1,0,N-1,l1,l2);
			}

		}

		ops_count--;
	}


	free(tree_data);
	free(data);
}