#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <math.h>
#include "linalg.h"

void scalar();
double (*generate_array())[3];
double calc_Distance(double (*locations)[3], int index1, int index2);
void writeToFile();
int readFile();

int main() {
    scalar();

    int n = 20;
    double (*LocArray)[3] = generate_array(n);
    int loc1 = 0;
    int loc2 = 2;
    double distance = calc_Distance(LocArray, loc1, loc2);
    printf("\nThe distance between the locations at index %d and %d is: %lf", loc1, loc2, distance);
    free(LocArray);


    writeToFile();
    readFile();
    return 0;
}

int readFile() {
   FILE *file;
   double array[5500];
   
   file = fopen("lotsOf.data", "rb");
   if (file == NULL) {
       perror("Error opening file");
       return 1;
   }

   size_t result = fread(array, sizeof(double), 5500, file);
   if (result != 5500) {
       perror("Error reading file");
   } else {
       for (int i = 0; i < 20; i++) {
           printf("Number %d: %f\n", i + 1, array[i]);
       }
   }
   
   fclose(file);
   return 0;
}

void writeToFile() {
    const gsl_rng_type * T;
    gsl_rng * r;
    double *longArray = malloc(sizeof(double)*5500);

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    int seed = 42;
    gsl_rng_set(r, seed); 

    int len = 5500;
    for (int i = 0; i < len; i++) {
        longArray[i] = gsl_rng_uniform(r);
    }

    FILE *f = fopen("lotsOf.data", "wb");
    int numOfSuccessfulElements = fwrite(longArray, sizeof(double), len, f);
    fclose(f);
    printf("\nNum of Successful Elements: %d \n", numOfSuccessfulElements);
    printf("Size of long Array: %ld\n", sizeof(double) * len);
    free(longArray);
}

double (*generate_array(int n))[3] {
    const gsl_rng_type * T;
    gsl_rng * r;
    double (*dynArray)[3] = malloc(sizeof(double[n][3]));

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    int seed = 42;
    gsl_rng_set(r, seed); 
    for (int i = 0; i < n; i++) {
        dynArray[i][0] = gsl_ran_gaussian(r, 1);
        dynArray[i][1] = gsl_ran_gaussian(r, 1); 
        dynArray[i][2] = gsl_ran_gaussian(r, 1);    
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < 3; ++j) {
            printf("dynArray[%d][%d] = %lf    ", i, j, dynArray[i][j]);
        }
        printf("\n");
    }
    gsl_rng_free(r);

    return dynArray;
}

double calc_Distance(double (*locations)[3], int index1, int index2) {
    double x_dist = locations[index1][0] - locations[index2][0];
    double y_dist = locations[index1][1] - locations[index2][1]; 
    double z_dist = locations[index1][2] - locations[index2][2];  
    return sqrt(x_dist * x_dist + y_dist * y_dist + z_dist * z_dist);
}

void scalar() {
    const gsl_rng_type * T;
    gsl_rng * r;
    int n;
    printf("The lenght of the vector should be: ");
    scanf("%d", &n);
    printf("\nThe lenght of the vector is: %d \n", n);

    double *arr1 = calloc(n, sizeof(double));
    double *arr2 = calloc(n, sizeof(double));

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    int seed = 42;
    gsl_rng_set(r, seed); 

    int i;
    for (i = 0; i < n; i++){
        arr1[i] = gsl_ran_gaussian(r, 1);
        arr2[i] = gsl_rng_uniform(r);
    }

    double scalar_product = scalarProduct(arr1, arr2, n);
    for (int j = 0; j < n; ++j) {
        if (j == 0) {
            printf("Vector 1: [%lf, ", arr1[j]);
        } else if (j == (n-1)) {
            printf(" %lf] \n", arr1[j]);
        } else {
            printf("%lf, ", arr1[j]);
        }
    }

    for (int k = 0; k < n; ++k) {
        if (k == 0) {
            printf("Vector 2: [%lf, ", arr2[k]);
        } else if (k == (n-1)) {
            printf(" %lf] \n", arr2[k]);
        } else {
            printf("%lf, ", arr2[k]);
        }
    }
    printf("\nThe scalar product is equal to: %lf \n \n", scalar_product);

    gsl_rng_free (r);
}