void bubblesort(unsigned long nel, int (*compare)(unsigned long i, unsigned long j), void (*swap)(unsigned long i, unsigned long j))
{
    int i = 0, ready = 0, lastswapright, lastswapleft, lborder, rborder;
    unsigned long left = 0, right = nel;
    right--;
    lborder = 0;
    rborder = nel - 1;
    do {
        ready = 0; 
       
        /* --> */
        for (i = left; i < right; i++) {
            if (compare(i + 1, i) < 0) { 
                swap(i + 1, i);
                ready = 1;           
                lastswapright = i;    
            }
        }
 

        int lsr = lastswapright;
        int lsl = lastswapleft;
        if ((ready == 0) && ((i == right) ||
            (lsl == left) || (rborder == 0) || (lborder == nel))) 
            goto exit;
 
        right = lastswapright; 
        ready = 0;
        lborder++;
 
 /* <-- */
        for (i = right; i > left; i--) { // <--
            if (compare(i - 1, i) > 0) {
                swap(i - 1, i);
                ready = 1;  
                lastswapleft = i;
            }
        }
 

        lsr = lastswapright;
        lsl = lastswapleft;
        if ((ready == 0) && ((lsr == right) ||
            (i == left) || (rborder == 0) || (lborder == nel)))
            goto exit;
 
        left = lastswapleft;
        rborder--;
    } while (ready);
    exit: ready = 0;
}