unsigned long binsearch(unsigned long nel, int (*compare)(unsigned long i))
{
    unsigned long low = 0, high, mid;
    high = nel-1;
    while (low < high-1) {
            mid = (low+high)/2;
            if (compare(mid)>0)
                high = mid;
            else
                if (compare(mid) < 0)
                    low = mid;
                else
                    if (compare(mid) == 0)
                        return mid;

    }
    if (compare(low) == 0) return low;
    else
        if (compare(high) == 0) return high;
        else
            return nel;
}