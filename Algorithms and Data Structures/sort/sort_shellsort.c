//http://sorting.at/
//


void shellsort(unsigned long nel, int (*compare)(unsigned long, unsigned long), void (*swap)(unsigned long, unsigned long))
{
    int i = 2, max_fib_index = 2;
    unsigned long f[nel];
    f[0] = 1;
    f[1] = 1;
    for (i; ;i++) {
        f[i] = f[i-1] + f[i-2];
        if (f[i] > nel) break;
        max_fib_index++; 
    }
    int gap = f[max_fib_index], ready = 0;

    int j;
    do {
           gap = f[max_fib_index];
	   i = 0;
	   ready = 0;

	   while ((i + gap) < nel){
                j = i + gap;
		while((j >= gap) && (compare(j, j - gap) < 0)) {
                        swap(j, j - gap);
                        j = j - gap;
		}
		i++;
	   }
	   max_fib_index--;
     } while (gap > 1);
}


//============================================================================


int compare(unsigned long i, unsigned long j) 
{ 
        if (i <= j) { 
                printf("COMPARE␣%d␣%d\n", i, j); 
        } else { 
                printf("COMPARE␣%d␣%d\n", j, i); 
        } 
 
        if (array[i] == array[j]) return 0; 
        return array[i] < array[j] ? -1 : 1; 
} 
 
void swap(unsigned long i, unsigned long j) 
{ 
        if (i <= j) { 
                printf("SWAP␣%d␣%d\n", i, j); 
        } else { 
                printf("SWAP␣%d␣%d\n", j, i); 
        } 
 
        int t = array[i]; 
        array[i] = array[j]; 
        array[j] = t; 
}