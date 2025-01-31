#include <stdio.h>

int
main()
{ 
    int m = 2;
    int *pm = &m;
    int **ppm = &pm;
    int ***pppm = &ppm;
    /* A pointer to a pointer */
    //int **ppm = &pm;

    printf("m + 1 = %i\n", m + 1);
    printf("(*pm) + 1 = %i\n", (*pm) + 1);
    printf("pm = %p\n", pm);
    printf("*(&m) + 1 = %i\n", *(&m) + 1);
    printf("---------------\n");
    printf("%i \n", pppm);
    printf("%i \n", *pppm);
    printf("%i \n", **pppm);
    printf("%i \n", ***pppm);
    return 0;
}
