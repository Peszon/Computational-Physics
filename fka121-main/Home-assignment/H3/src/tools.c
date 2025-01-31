#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <complex.h>

// #include "/home/felixpersson/kompfys/fka121-main/Exercies/include/mechanics_utilities.h"

void
remove_element_at_index(
        double *list,
        int index,
        int len 
        )
{
    list[index] = list[len - 1];
    list[len - 1] = NAN;  
}

void
remove_list_at_index(
        double **list,
        int index,
        int n,
        int m 
        )
{
    free(list[index]);
    list[index] = list[n - 1];
    list[n - 1] = malloc(sizeof(double) * m);
}

void
elementwise_addition(
                     double *res,
                     double *v1,
                     double *v2,
                     unsigned int len
                    )
{
    for(int i = 0; i < len; ++i){
	    res[i] = v1[i] + v2[i];
    }
}

void
elementwise_multiplication(
                           double *res,
                           double *v1,
                           double *v2,
                           unsigned int len
                          )
{
    for(int i = 0; i < len; ++i){
	    res[i] = v1[i] * v2[i];
    }
}

void
addition_with_constant(
                       double *res,
                       double *v,
                       double constant,
                       unsigned int len)
{
    for(int i = 0; i < len; ++i){
	    res[i] = v[i] + constant;
    }
}

void
multiplication_with_constant(
                             double *res,
                             double *v,
                             double constant,
                             unsigned int len)
{
    for(int i = 0; i < len; ++i){
	    res[i] = v[i] * constant;
    }
}

double
dot_product(
            double *v1,
            double *v2,
            unsigned int len
           )
{
    double res = 0;
    for(int i = 0; i < len; ++i){
	    res += v1[i] * v2[i];
    }
    return res;
}

double *
linspace(
    double start,
    double stop,
    int N)
{
    double *array = malloc(N * sizeof(double));

    double delta = (stop - start) / N;

    for (int i = 0; i < N; i++) {
        array[i] = start + i * delta;
    }

    return array;
}

double **
create_2D_array(
                unsigned int row_size,
                unsigned int column_size
               )
{
    double **array = malloc(row_size * sizeof(double *));
    if (array == NULL) {
        perror("Failed to allocate memory");
        return NULL;
    }

    for (int i = 0; i < row_size; i++) {
        array[i] = malloc(column_size * sizeof(double));
        if (array[i] == NULL) {
            perror("Failed to allocate memory");

            for (unsigned int j = 0; j < i; j++) {
                free(array[j]);
            }

            free(array);
            return NULL;
        }
    }
    return array;
}

double ***create_3D_array(
    unsigned int row_size,
    unsigned int column_size,
    unsigned int depth_size)
{
    // Allocate memory for the outermost array (depth level)
    double ***array = malloc(depth_size * sizeof(double **));
    if (array == NULL) {
        perror("Failed to allocate memory for depth level");
        return NULL;
    }

    for (unsigned int i = 0; i < depth_size; i++) {
        // Allocate memory for each row-level array
        array[i] = malloc(row_size * sizeof(double *));
        if (array[i] == NULL) {
            perror("Failed to allocate memory for row level");

            // Free already allocated memory
            for (unsigned int j = 0; j < i; j++) {
                free(array[j]);
            }
            free(array);
            return NULL;
        }

        for (unsigned int j = 0; j < row_size; j++) {
            // Allocate memory for each column
            array[i][j] = malloc(column_size * sizeof(double));
            if (array[i][j] == NULL) {
                perror("Failed to allocate memory for column level");

                // Free already allocated memory
                for (unsigned int k = 0; k < j; k++) {
                    free(array[i][k]);
                }
                free(array[i]);

                for (unsigned int k = 0; k < i; k++) {
                    for (unsigned int l = 0; l < row_size; l++) {
                        free(array[k][l]);
                    }
                    free(array[k]);
                }
                free(array);
                return NULL;
            }
        }
    }

    return array;
}

void
destroy_2D_array(
                 double **array,
                 unsigned int n
                )

{
    for (unsigned int i = 0; i < n; i++) {
        free(array[i]);
    }
    
    free(array);
}

void 
destroy_3D_array(
                double ***array, 
                unsigned int depth_size, 
                unsigned int row_size) 
{
    for (unsigned int i = 0; i < depth_size; i++) {
        for (unsigned int j = 0; j < row_size; j++) {
            free(array[i][j]);
        }
        free(array[i]);
    }
    free(array);
}

void
matrix_vector_multiplication(
                             double *result,
                             double **A,
                             double *b,
                             unsigned int n,
                             unsigned int m
                            )
{
    for (int i = 0; i < n; i++) {
        result[i] = dot_product(A[i], b, m);
    }
}

void
matrix_matrix_multiplication(
                             double **result,
                             double **A,
                             double **B,
                             unsigned int n,
                             unsigned int m,
                             unsigned int k
                            )
{
    for (int i = 0; i < n; i++) {
        for (int l = 0; l < k; l++) {
            result[i][l] = 0;
            for (int j = 0; j < m; j++) {
                result[i][l] += A[i][j] * B[j][l];
            }
        }
    }
}

double
vector_norm(
            double *v1,
            unsigned int len
           )
{
    double res = 0; 
    for (int i = 0; i < len; i++) {
        res += pow(v1[i], 2);
    }
    return sqrt(res);    
}


void
normalize_vector(
                 double *v1,
                 unsigned int len
                )
{
    double norm = vector_norm(v1, len);
    double norm_factor = 1 / norm;
    multiplication_with_constant(v1, v1, norm_factor, len);
}

double
average(
        double *v1,
        unsigned int len
       )
{
    double sum = 0;
    for (int i = 0; i < len; i++) {
        sum += v1[i];
    }
    return sum / len;
}


double
standard_deviation(
                       double *v1,
                       unsigned int len
                  )
{
    double mu = average(v1, len);
    double residue = 0;
    for (int i = 0; i < len; i++) {
        residue += pow(v1[i] - mu, 2);
    }
    return sqrt(residue/len);
}

double
variance(
                       double *v1,
                       unsigned int len
                  )
{
    double mu = average(v1, len);
    double residue = 0;
    for (int i = 0; i < len; i++) {
        residue += pow(v1[i] - mu, 2);
    }
    return residue/len;
}

double
distance_between_vectors(
                         double *v1,
                         double *v2,
                         unsigned int len
                        )
{
    double distance[len];
    multiplication_with_constant(v1, v1, -1, len);
    elementwise_addition(distance, v1, v2, len);
    return vector_norm(distance, len);
}

void
cumulative_integration(
                       double *res,
                       double *v,
                       double dx,
                       unsigned int v_len
                      )
{
    res[0] = 0;
    for (int i = 1; i < v_len; i++) {
        res[i] = res[i-1] + dx * (v[i-1] + v[i]) / 2; 
    }
}


void
write_xyz(
          FILE *fp,
          char *symbol,
          double **positions,
          double **velocities,
          double alat,
          int natoms)
{
    fprintf(fp, "%i\nLattice=\"%f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f\" ", natoms, alat, alat, alat);
    fprintf(fp, "Properties=species:S:1:pos:R:3:vel:R:3 pbc=\"T T T\"\n");
    for(int i = 0; i < natoms; ++i){
        fprintf(fp, "%s %lf %lf %lf %lf %lf %lf\n",
                                                symbol, 
                                                positions[i][0], 
                                                positions[i][1], 
                                                positions[i][2], 
                                                velocities[i][0],
                                                velocities[i][1],
                                                velocities[i][2]);
    }
}

void fft_freq(
          double *res,
              int n,
              double timestep)
{
    if (n % 2 == 0) {
        for (int i = 0; i < (n/2); i++) {
            res[i] = i;
        }
        for (int i = 0; i < n/2; i++) {
            res[i + n/2] = -n/2 + i;
        }
    } else {
        for (int i = 0; i <= ((n-1)/2); i++) {
            res[i] = i;
        }
        for (int i = 0; i < (n-1)/2; i++) {
            res[i + (n+1)/2] = -(n-1)/2 + i;
        }
    } 
    double factor = 2 * M_PI / (n * timestep);
    multiplication_with_constant(res, res, factor, n);
}

/* Freely given functions */
void
skip_line(FILE *fp)
{
    int c;
    while (c = fgetc(fp), c != '\n' && c != EOF);
}

void
read_xyz(
         FILE *fp,
         char *symbol,
         double **positions,
         double **velocities,
         double *alat)
{
    int natoms;
    if(fscanf(fp, "%i\nLattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 %lf\" ", &natoms, alat, alat, alat) == 0){
        perror("Error");
    }
    skip_line(fp);
    for(int i = 0; i < natoms; ++i){
        int successful1 = fscanf(fp, "%s %lf %lf %lf ",
                symbol, &positions[i][0], &positions[i][1], &positions[i][2]);
        int successful2 = fscanf(fp, "%lf %lf %lf\n",
                &velocities[i][0], &velocities[i][1], &velocities[i][2]);
        if ((successful1 == 0) || (successful2 == 0)) {
            perror("Failed to read extxyz file.");
        }
    }
}

void powerspectrum(
           double *res,
           double *signal,
           int n,
           double timestep)
{
    /* Declaration of variables */
    double *complex_coefficient = malloc(sizeof(double) * 2*n); // array for the complex fft data
    double *data_cp = malloc(sizeof(double) * n);

    /*make copy of data to avoid messing with data in the transform*/
    for (int i = 0; i < n; i++) {
    data_cp[i] = signal[i];
    }

    /* Declare wavetable and workspace for fft */
    gsl_fft_real_wavetable *real;
    gsl_fft_real_workspace *work;

    /* Allocate space for wavetable and workspace for fft */
    work = gsl_fft_real_workspace_alloc(n);
    real = gsl_fft_real_wavetable_alloc(n);

    /* Do the fft*/
    gsl_fft_real_transform(data_cp, 1, n, real, work);

    /* Unpack the output into array with alternating real and imaginary part */
    gsl_fft_halfcomplex_unpack(data_cp, complex_coefficient,1,n);

    /*fill the output powspec_data with the powerspectrum */
    for (int i = 0; i < n; i++) {
    res[i] = (complex_coefficient[2*i]*complex_coefficient[2*i]+complex_coefficient[2*i+1]*complex_coefficient[2*i+1]);
    res[i] *= timestep / n;
    }

    /* Free memory of wavetable and workspace */
    gsl_fft_real_wavetable_free(real);
    gsl_fft_real_workspace_free(work);
    free(complex_coefficient);
    free(data_cp);
}

void
write_2D_array_to_csv(
    const char *filename, 
    double **array, 
    int rows, 
    int cols) 
{
    FILE *file = fopen(filename, "w");
    if (file ==  NULL) {
        perror("Error: invalid filename/path");
        return;
    }
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            fprintf(file, "%f", array[i][j]);
            if (j != cols-1) {
                fprintf(file, ",");
            }
        }
        fprintf(file, "\n");
    }
    fclose(file);
    return;
}

void
write_1D_array_to_csv(
    const char *filename, 
    double *array, 
    int len)
{
    FILE *file = fopen(filename, "w");
    if (file ==  NULL) {
        perror("Error: invalid filename/path");
        return;
    }
    for (int i = 0; i < len; ++i) {
        fprintf(file, "%f", array[i]);
        if (i != len-1) {
            fprintf(file, ",");
        }
        fprintf(file, "\n");
    }
    fclose(file);
    return;
}

double autocorrelation(
            double *data,
            int data_len,
            int time_lag_ind)
{
    double *normalized_data = malloc(data_len * sizeof(double));
    double data_mean = average(data, data_len);

    addition_with_constant(normalized_data, data, -data_mean, data_len);

    double sum = 0.0, sum_sq = 0.0, sum_prod = 0.0;
    int valid_pairs = data_len - time_lag_ind;

    for (int i = 0; i < data_len; i++) {
        sum += normalized_data[i];
        sum_sq += normalized_data[i] * normalized_data[i];
        if (i + time_lag_ind < data_len) {
            sum_prod += normalized_data[i] * normalized_data[i + time_lag_ind];
        }
    }

    free(normalized_data);

    double mean = sum / data_len;
    double mean_sq = sum_sq / data_len;
    double mean_prod = sum_prod / valid_pairs;

    return (mean_prod - (mean * mean)) / (mean_sq - (mean * mean));
}

double block_average(double *data,
                int data_len,
                int block_size
                )
    {
    double *normalized_data = malloc(data_len * sizeof(double));
    double data_mean = average(data, data_len);

    addition_with_constant(normalized_data, data, -data_mean, data_len);

    int num_blocks = (int) data_len / block_size;
    double *block_average_arr = malloc(num_blocks * sizeof(double));

    for (int i = 0; i < num_blocks; i++) {
        double block_sum = 0;
        for (int j = 0; j < block_size; j++) {
            block_sum += normalized_data[j + block_size * i];
        }
        block_average_arr[i] = block_sum / block_size; 
    }

    double std_block_average = standard_deviation(block_average_arr, num_blocks);
    double std_data = standard_deviation(normalized_data, data_len);

    free(block_average_arr);
    free(normalized_data);

    return block_size * pow(std_block_average,2) / pow(std_data,2);
}
