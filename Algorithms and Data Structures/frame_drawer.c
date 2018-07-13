#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void top_bot_frame (int width)
{
    int i;
    for (i = 0; i < width; i++) printf("*");
    printf("\n");
}

void free_space (int height, int width, char *str)
{
    int i, j, l = height/2 - 1;
    for (i = 0; i < l; i++) {
        printf("*");
        for (j = 0; j < width - 2; j++) printf(" ");
        printf("*");
        printf("\n");
    }
}

int main(int argc, char **argv)
{
    int h, w, k, i;
    char string[100];

    //h = argv[1];
    //w = argv[2];
    scanf("%d", &h);
    scanf("%d", &w);
    scanf("%s", string);

    /*if ( < 4) {
            printf ("Usage: frame <height> <width> <text>");
            return -1;
    }*/

    //char *string = argv[3];

    int len = strlen(string);

    if (w < len + 2) {
        printf("Error\n");
        return -1;
    }

    k = (w - len) / 2 - 1;

    top_bot_frame (w);
    free_space (h, w, string);

    printf("*");
    for (i = 0; i < k; i++) printf(" ");
    printf("%s", string);
    for (i = 0; i < k; i++) printf(" ");
    if ((w % 2 == 0) && (len % 2 != 0)) printf(" ");
    printf("*");
    printf("\n");

    free_space(h, w, string);
    top_bot_frame (w);

    return 0;
}