#include "linalg.h"

double scalarProduct(double *arr1, double *arr2, int len) {
    double scalar_product = 0;
    for (int i = 0; i < len; i++) {
        scalar_product += arr1[i] * arr2[i];
    }

    return scalar_product;
}