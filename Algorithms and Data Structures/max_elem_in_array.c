int maxarray (void *base, unsigned long nel, unsigned long width,
              int (*compare) (void *a, void *b))
{
    int i = 0, j = 0;
    void  *maxel = base;
    i = 0;
    j = 0;
    for (i = 0; i < nel; i++) {
        if (compare((base+width*i), maxel) > 0) {
            maxel = base+width*i;
            j = i;
        }
    }
    return j;
}