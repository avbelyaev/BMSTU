unsigned long peak(unsigned long nel, int (*less)(unsigned long i, unsigned long j))
{
    unsigned long low, high;
    low = 0;
    high = nel;

    while (low <= high) {
        unsigned long mid = (high - low) / 2 + low;
        if (less(mid, mid-1)) high = mid;
        else
            if (less(mid, mid+1)) low = mid;
        else
            return mid;
        }
}