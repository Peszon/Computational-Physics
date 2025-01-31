#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "tools.h"

#define PI 3.141592653589


typedef struct {
    double integral;
    double error;
} result_t;

typedef struct{
    double weight;
    double function_value;
    int accepted;
} result_d;


/* **********************************************
 *
 * Perform Monte Carlo integration of the given
 * integral without using importance sampling.
 *
 * Parameters
 * ----------
 *  N - Number of points to sample
 *  k - GSL random number generator object
 *
 * Returns
 * -------
 *  Struct with the result
 *
 * **********************************************/
result_t MC_without_importance_sampling(int N, gsl_rng *k)
{
    result_t result;

    double sum = 0;
    double sum2 = 0;
    for(int i = 0; i < N; i++) {
        double x_i = gsl_rng_uniform(k);
        sum += x_i * (1 - x_i);
        sum2 += pow(x_i*(1-x_i), 2); 
    }
    result.integral = sum/N;
    result.error = pow(sum2/N - pow(sum/N, 2), 1./2.) / pow(N, 1./2.);
    return result;
}

/* **********************************************
 *
 * Perform Monte Carlo integration of the given
 * integral using importance sampling.
 *
 * Parameters
 * ----------
 *  N - Number of points to sample
 *  k - GSL random number generator object
 *
 * Returns
 * -------
 *  Struct with the result
 *
 * **********************************************/
result_t MC_with_importance_sampling(int N, gsl_rng *k)
{
    result_t result;
    double sum = 0;
    double sum2 = 0;
    for(int i = 0; i < N; i++) {
        double chi_i = gsl_rng_uniform(k);
        double x_i = 1/PI * acos(1 - 2*chi_i);

        sum += x_i * (1 - x_i) / (sin(x_i*PI) * PI / 2);
        sum2 += pow(x_i*(1-x_i) / (sin(x_i*PI)* PI / 2), 2); 
    }
    result.integral = sum/N;
    result.error = pow(sum2/N - pow(sum/N, 2), 1./2.) / pow(N, 1./2.);
    return result;
}

void task3() {
        int N[] = {10, 100, 1000, 10000};
    result_t results_with;
    result_t results_without;

    const gsl_rng_type * T;
    gsl_rng * r;

    int i, n = 10;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    int seed = 42;
    gsl_rng_set(r, seed); 

    for (int i = 0; i < 4; i++) {
        result_t results_without = MC_without_importance_sampling(N[i], r);
        result_t results_with = MC_with_importance_sampling(N[i], r);
        printf("N = %d, Integral with importance sampling: %lf +- %lf, Integral without importance sampling: %lf +- %lf \n", N[i], results_with.integral, results_with.error, results_without.integral, results_without.error);
    }
    printf("Correct answer: %lf \n", 1./6.);

    gsl_rng_free (r);
}


double weight(
    double *x)
{
    return exp(-1.*(pow(x[0],2.) + pow(x[1],2.) + pow(x[2],2.)));
}

double function(
                double *x)
{

    return pow(x[0], 2) + pow(x[0]*x[1], 2) + pow(x[0]*x[1]*x[2], 2);
}

 result_d MCMC_step_displace_all(
                                 double *x,
                                 double delta,
                                 gsl_rng *k)
{
    result_d result;
    double x_new[3];
    x_new[0] = x[0] + delta * (gsl_rng_uniform(k) - 0.5);
    x_new[1] = x[1] + delta * (gsl_rng_uniform(k) - 0.5);
    x_new[2] = x[2] + delta * (gsl_rng_uniform(k) - 0.5);

    double accept_value = gsl_rng_uniform(k);
    if (weight(x_new)/weight(x) > accept_value) {
        result.weight = weight(x_new);
        result.function_value = function(x_new);
        result.accepted = 1;

        x[0] = x_new[0];
        x[1] = x_new[1];
        x[2] = x_new[2];
    } else {
        result.weight = weight(x);
        result.function_value = function(x);
        result.accepted = 0;
    }

    return result;
}

void task5() {
    int samples = 10000000;
    int accepted = 0;
    double x[3] = {0};
    double delta = 2;
    double sum = 0;

    const gsl_rng_type * T;
    gsl_rng * k;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    k = gsl_rng_alloc(T);

    int seed = 42;
    gsl_rng_set(k, seed); 

    for (int i = 0; i<samples; i++) {
        result_d result = MCMC_step_displace_all(x, delta, k);
        sum += result.function_value;
        accepted += result.accepted;
    }

    printf("The integral was evaluated to: %lf \n", sum/samples);
    printf("The acceptance ratio of the proposed step is: %lf \n", (double) accepted/samples);
    gsl_rng_free(k);
}

void task7() {
    // Opening the file and reading the data.
    FILE *file;
    char buffer[256]; // Buffer to store each line
    
    file = fopen("MC.txt", "r");
    if (file == NULL) {
        perror("Error opening file");
    }

    int data_length = 1000000; 
    double data[data_length];
    int i = 0;
    while (fgets(buffer, sizeof(buffer), file) != NULL) {
        data[i] = strtod(buffer,NULL);
        i++; 
    }

    // Calculate the statistical inefficiency using the block average method. 
    int block_length_min = 1;
    int block_length_max = 800;
    double statistical_inefficiency_block[block_length_max-block_length_min];
    for (int i = 0; i < block_length_max-block_length_min; i++) {
        if (i % 50 == 0) {
            printf("i: %d, block_stat_ineff: %lf \n", i, block_average(data, data_length, (i+block_length_min)));
        }
        statistical_inefficiency_block[i] = block_average(data, data_length, (i+block_length_min)*5);
    }
    write_1D_array_to_csv("statistical_ineff_block_data.dat", statistical_inefficiency_block, block_length_max-block_length_min);

    // Calculate the statistical inefficiency using the autocorrelation function.
    for (int k = 1; k < 2000; k++) {
        if (autocorrelation(data, data_length, k) < 0.1) {
            printf("statistical ineff using autocorrelation func: %d", k);
            break;
        }
        if (k % 10 == 0) {
            printf("k: %d, autocorrelation_stat_ineff: %lf \n \n", k, autocorrelation(data, data_length, k));
        }
    }

    fclose(file);
}

int main() {
    task3();
    task5();
    return 1;
}
