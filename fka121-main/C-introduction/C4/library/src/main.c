#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "/home/felixpersson/kompfys/fka121-main/C-introduction/C4/library/include/tools.h"

int main()
{
    /* *******************************
     * 
     * Here is some code that is
     * using vector.c
     *
     * For an example look at the
     * third exercise
     *
     * ******************************/
    double **A = create_2D_array(5,5);
    double **B = create_2D_array(5,5);
    A[4][4] = 5;
    B[4][4] = 10;

    double **C = create_2D_array(5,5);
    matrix_matrix_multiplication(C, A, B, 5,5,5);
    printf("C[0][0] = %lf\n", C[4][4]);
    return 0;
}
