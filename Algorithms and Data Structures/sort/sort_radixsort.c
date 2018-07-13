#include <stdio.h>
#include <stdlib.h>
#define max_day 31
#define min_day 1
#define max_month 12
#define min_month 1
#define max_year 2030
#define min_year 1970

struct Date { int day, month, year; };

int main()
{
    int i, n;
    scanf ("%d", &n);

    struct Date arr[n];
    struct Date temp_arr[n];

    void radixsort (int marker, int range)
    {
        int i, temp, bucket[100] = { 0 };
        if (marker == 1) {
            for (i = 0; i < n; i++) {
                temp = arr[i].day - min_day;
                bucket[temp]++;
            }
        }
        else if (marker == 2) {
                for (i = 0; i < n; i++) {
                    temp = arr[i].month - min_month;
                    bucket[temp]++;
                }
            }
            else {
                for (i = 0; i < n; i++) {
                temp = arr[i].year - min_year;
                    bucket[temp]++;
                }
            }
        for (i = 1; i < range; i++) bucket[i] += bucket[i - 1];

        if (marker == 1) {
            for (i = n - 1; i >= 0; i--) {
                temp = arr[i].day - min_day;
                temp_arr[--bucket[temp]] = arr[i];
            }
        }
        else
            if (marker == 2) {
                for (i = n - 1; i >= 0; i--) {
                    temp = arr[i].month - min_month;
                    temp_arr[--bucket[temp]] = arr[i];
                }
            }
            else {
                for (i = n - 1; i >= 0; i--) {
                    temp = arr[i].year - min_year;
                    temp_arr[--bucket[temp]] = arr[i];
                }
            }

        for (i = 0; i < n; i++) arr[i] = temp_arr[i];
    }

    for (i = 0; i < n; i++) {
        scanf ("%d", &arr[i].year);
        scanf ("%d", &arr[i].month);
        scanf ("%d", &arr[i].day);
    }

    int band = max_day;/
    int pointer = 1;
    radixsort (pointer, band);

    band = max_month;
    pointer = 2;
    radixsort (pointer, band);

    band = max_year - min_year + 1;
    pointer = 3;
    radixsort (pointer, band);

    for (i = 0; i < n; i++) {
        printf ("%d ", arr[i].year);

            printf ("%d ", arr[i].month);

            printf ("%d\n", arr[i].day);
    }
    return 0;
}