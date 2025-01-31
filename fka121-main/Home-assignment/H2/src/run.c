#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>

#include "tools.h"
#include "lattice_utils.h"

void task1();
double newton_raphson(int maxmitr, double allerr, double temperature, double delta_E);

void task2();

void task3();

void task3_optimized();
double calculate_stat_ineff_autocorr(double *data, int data_length);
double calculate_stat_ineff_block(double *data, int data_length);

int
run(
    int argc,
    char *argv[]
   )
{
    task3_optimized();
    return 0;
}

void task1(){
    double E_CC = -0.436;  // eV
    double E_ZZ = -0.113;  // eV
    double E_CZ = -0.294;  // eV

    double delta_E = E_CC + E_ZZ - 2 * E_CZ;
    int samples = 100000;
    double *temperature = linspace(1,1001, samples);    // K

    double results_P[samples];
    double results_U[samples];

    for (int i = 0; i < samples; i++) {
        results_P[i] = newton_raphson(2000, 1e-9, temperature[i], delta_E);
    }

    write_1D_array_to_csv("P_of_T.dat", results_P, samples);
    for (int i = 0; i < samples; i++) {
        results_U[i] = -pow(results_P[i], 2);
    }
    write_1D_array_to_csv("U_of_T.dat", results_U, samples);
}

float f(double P, double T, double delta_E)
{   
    return k_B * T * log((1+P)/(1-P)) - 4*delta_E*P;
}

float df (double P, double T, double delta_E)
{
    return k_B * T*2 / (1 - pow(P,2)) - 4 * delta_E;
}


double newton_raphson(int maxmitr, double allerr, double temperature, double delta_E)
{
    double h;
    double P0 = 0.9999999999;
    double P1;
    for (int itr=1; itr<=maxmitr; itr++)
    {
        if (df(P0, temperature, delta_E) < 1e-8) {
            //printf("After %3d iterations, root = %8.6f, temperature: %lf\n", maxmitr, P0, temperature - 273.15);
            break;
        } else {
            h = f(P0, temperature, delta_E)/df(P0, temperature, delta_E);
        }
        P1 = P0 - h; 
        if (P1 > 1) {
            P1 = 0.9999999;
        }
        if ((fabs(h) < allerr))
        {
            return P1;
        }
        P0 = P1;
    }
    return 1e100;
}

void
task2() {
    // Setting up GSL rng.
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    int seed = 42;
    gsl_rng_set(r, seed); 

    // Creating the lattices.
    int lattice_size = 10;
    double ***lattice_A = create_cold_lattice(lattice_size, 1.);
    double ***lattice_B = create_cold_lattice(lattice_size, 2.);

    // Monte carlo simulation settings and data.
    int eq_steps = 1e6;
    //int meas_steps = 1000;
    double temp = 400;     // K
    double E_CC = -0.436;  // eV
    double E_ZZ = -0.113;  // eV
    double E_CZ = -0.294;  // eV
    double energy = calculate_BCC_lattice_energy(lattice_A, lattice_B, lattice_size, E_ZZ, E_CC, E_CZ);
    double *energy_history = malloc(sizeof(double)*eq_steps); 
    double *long_range_order_history = malloc(sizeof(double)*eq_steps); 
    double *short_range_order_history = malloc(sizeof(double)*eq_steps); 

    for (int i = 0; i < eq_steps; i++) {
        int outcome = metropolis_one_step(lattice_A, lattice_B, lattice_size, E_ZZ, E_CC, E_CZ, temp, &energy, r);
        energy_history[i] = energy;
        long_range_order_history[i] = calculate_long_range_order(lattice_A, lattice_size);
        short_range_order_history[i] = calculate_short_range_order(lattice_A, lattice_B, lattice_size);
        if (i % 10000 == 0) {
            printf("iteration %d: Outcome = %d, Energy = %lf \n", i, outcome, energy - energy_history[0]);
        }
    }

    write_1D_array_to_csv("E_hist_400K.dat", energy_history, eq_steps);
    write_1D_array_to_csv("long_hist_1000K.dat", long_range_order_history, eq_steps);
    write_1D_array_to_csv("short_hist_1000K.dat", short_range_order_history, eq_steps);

    free(energy_history);
    free(long_range_order_history);
    free(short_range_order_history);
    destroy_3D_array(lattice_A, lattice_size, lattice_size);
    destroy_3D_array(lattice_B, lattice_size, lattice_size);
}

void
task3() {
    // Setting up GSL rng.
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    int seed = 42;
    gsl_rng_set(r, seed); 

    // Creating the lattices.
    int lattice_size = 10;
    double ***lattice_A = create_cold_lattice(lattice_size, 1.);
    double ***lattice_B = create_cold_lattice(lattice_size, 2.);

    // Monte carlo sampling settings.
    int eq_steps = 5e4;
    int meas_steps = 2e5;
    double E_CC = -0.436;  // eV
    double E_ZZ = -0.113;  // eV
    double E_CZ = -0.294;  // eV
    double energy = calculate_BCC_lattice_energy(lattice_A, lattice_B, lattice_size, E_ZZ, E_CC, E_CZ);
    int temperature_steps = 3500;
    double *temperature = linspace(300, 1000, temperature_steps);
    int outcome;

    // Data saving arrays
    double **energy_history = create_2D_array(temperature_steps, meas_steps);
    double **long_range_order_history = create_2D_array(temperature_steps, meas_steps);
    double **short_range_order_history = create_2D_array(temperature_steps, meas_steps);


    // Sampling proccedure
    for (int i = 0; i < temperature_steps; i++) {

        // Equilibration samples
        for (int j = 0; j < eq_steps; j++) {
            outcome = metropolis_one_step(lattice_A, lattice_B, lattice_size, E_ZZ, E_CC, E_CZ, temperature[i], &energy, r);
            
            //Progress tracker
            if (j % 10000 == 0) {
                printf("j: %d, outcome: %d \n", j, outcome);
            }
        }

        // Measurment samples
        for (int j = 0; j < meas_steps; j++) {
            outcome = metropolis_one_step(lattice_A, lattice_B, lattice_size, E_ZZ, E_CC, E_CZ, temperature[i], &energy, r);
            energy_history[i][j] = energy;
            long_range_order_history[i][j] = calculate_long_range_order(lattice_A, lattice_size);
            short_range_order_history[i][j] = calculate_short_range_order(lattice_A, lattice_B, lattice_size);

            //Progress tracker
            if (j % 10000 == 0) {
                printf("j: %d, outcome: %d \n", j, outcome);
            }
        }

        //Progress tracker
        if (i % 10 == 0) {
            printf("i: %d, outcome: %d \n", i, outcome);
        }
    }

    // Saving data and freeing memory
    write_2D_array_to_csv("U_hist_300_1000K.dat", energy_history, temperature_steps, meas_steps);
    write_2D_array_to_csv("long_hist_300_1000K.dat", long_range_order_history, temperature_steps, meas_steps);
    write_2D_array_to_csv("short_hist_300_1000K.dat", short_range_order_history, temperature_steps, meas_steps);

    destroy_2D_array(energy_history, temperature_steps);
    destroy_2D_array(long_range_order_history, temperature_steps);
    destroy_2D_array(short_range_order_history, temperature_steps);
    destroy_3D_array(lattice_A, lattice_size, lattice_size);
    destroy_3D_array(lattice_B, lattice_size, lattice_size);
}

void task3_optimized() {
    // Simulation parameters
    int lattice_size = 10;
    double E_CC = -0.436;  // eV
    double E_ZZ = -0.113;  // eV
    double E_CZ = -0.294;  // eV

    // Temperature settings
    int temperature_steps = 350;
    double *temperature = linspace(300.0, 1000.0, temperature_steps); // K

    // Monte Carlo settings
    int eq_steps = 200000;
    int meas_steps = 300000;

    // Initialize the arrays used for saving
    double *energy_history = malloc(temperature_steps * sizeof(double));
    double *energy_corr_history = malloc(temperature_steps * sizeof(double));
    double *energy_block_history = malloc(temperature_steps * sizeof(double));

    double *capacity_history = malloc(temperature_steps * sizeof(double));
    double *capacity_corr_history = malloc(temperature_steps * sizeof(double));
    double *capacity_block_history = malloc(temperature_steps * sizeof(double));

    double *long_order_history = malloc(temperature_steps * sizeof(double));
    double *long_order_corr_history = malloc(temperature_steps * sizeof(double));
    double *long_order_block_history = malloc(temperature_steps * sizeof(double));

    double *short_order_history = malloc(temperature_steps * sizeof(double));
    double *short_order_corr_history = malloc(temperature_steps * sizeof(double));
    double *short_order_block_history = malloc(temperature_steps * sizeof(double));

    // Check for allocation success for each array (optional, but recommended)
    if (!energy_history || !energy_corr_history || !energy_block_history ||
        !capacity_history || !capacity_corr_history || !capacity_block_history ||
        !long_order_history || !long_order_corr_history || !long_order_block_history ||
        !short_order_history || !short_order_corr_history || !short_order_block_history) {
        fprintf(stderr, "Memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }

    // Parallelize over temperature steps
    #pragma omp parallel
    {
        // Initialize thread-local RNG
        gsl_rng_env_setup();
        const gsl_rng_type * T = gsl_rng_default;
        gsl_rng *r = gsl_rng_alloc(T);
        unsigned long seed = 42 + omp_get_thread_num();
        gsl_rng_set(r, seed);

        // Buffers to store measurement data for a single temperature
        double *energy_buffer = malloc(sizeof(double) * meas_steps);
        double *long_order_buffer = malloc(sizeof(double) * meas_steps);
        double *short_order_buffer = malloc(sizeof(double) * meas_steps);

        // Ensure memory allocation was successful
        if (!energy_buffer || !long_order_buffer || !short_order_buffer) {
            fprintf(stderr, "Memory allocation failed in thread %d.\n", omp_get_thread_num());
            exit(EXIT_FAILURE);
        }

        // Iterate over temperature steps
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < temperature_steps; i++) {
            double temp = temperature[i];

            double ***lattice_A = create_cold_lattice(lattice_size, 1.0);
            double ***lattice_B = create_cold_lattice(lattice_size, 2.0);

            double energy = calculate_BCC_lattice_energy(lattice_A, lattice_B, lattice_size, E_ZZ, E_CC, E_CZ);

            // Equilibration steps
            for (int j = 0; j < eq_steps; j++) {
                metropolis_one_step_optimized(lattice_A, lattice_B, lattice_size, E_CC, E_ZZ, E_CZ, temp, &energy, r);
                if (j % 50000 == 0 && omp_get_thread_num() == 0) {
                    printf("Temperature %d/%d (%lf K): Equilibration step %d\n", i+1, temperature_steps, temp, j);
                }
            }

            // Measurement steps
            for (int j = 0; j < meas_steps; j++) {
                metropolis_one_step_optimized(lattice_A, lattice_B, lattice_size, E_CC, E_ZZ, E_CZ, temp, &energy, r);
                energy_buffer[j] = energy;
                long_order_buffer[j] = calculate_long_range_order(lattice_A, lattice_size);
                short_order_buffer[j] = calculate_short_range_order(lattice_A, lattice_B, lattice_size);

                if (j % 50000 == 0 && omp_get_thread_num() == 0) {
                    printf("Temperature %d/%d (%lf K, %lf eV): Measurement step %d\n", i+1, temperature_steps, temp, energy, j);
                }
            }

            // Save measurement average data to array (synchronized)
            energy_history[i] = average(energy_buffer, meas_steps);
            capacity_history[i] = 1. / (k_B * pow(temp, 2)) * variance(energy_buffer, meas_steps);
            long_order_history[i] = average(long_order_buffer, meas_steps);
            short_order_history[i] = average(short_order_buffer, meas_steps);

            if (omp_get_thread_num() == 0) {
                printf("Completed %d/%d temperature steps. Average energy: %lf eV \n", i+1, temperature_steps, energy_history[i]);
            }

            energy_corr_history[i] = calculate_stat_ineff_autocorr(energy_buffer, meas_steps);
            long_order_corr_history[i] = calculate_stat_ineff_autocorr(long_order_buffer, meas_steps);
            short_order_corr_history[i] = calculate_stat_ineff_autocorr(short_order_buffer, meas_steps);

            if (omp_get_thread_num() == 0) {
                printf("Completed %d/%d temperature steps. Average energy: %lf eV \n", i+1, temperature_steps, energy_history[i]);
            }

            energy_block_history[i] = calculate_stat_ineff_block(energy_buffer, meas_steps);
            long_order_block_history[i] = calculate_stat_ineff_block(long_order_buffer, meas_steps);
            short_order_block_history[i] = calculate_stat_ineff_block(short_order_buffer, meas_steps);
            

            // Optional: Progress tracker
            if (omp_get_thread_num() == 0) {
                printf("Completed %d/%d temperature steps. Average energy: %lf eV \n", i+1, temperature_steps, energy_history[i]);
            }

            destroy_3D_array(lattice_A, lattice_size, lattice_size);
            destroy_3D_array(lattice_B, lattice_size, lattice_size);
        }

        // Free thread-local resources
        free(energy_buffer);
        free(long_order_buffer);
        free(short_order_buffer);
        gsl_rng_free(r);
    } // End of parallel region

    //Write to 
    write_1D_array_to_csv("U_hist_300_1000K.dat", energy_history, temperature_steps);
    write_1D_array_to_csv("U_corr_hist_300_1000K.dat", energy_corr_history, temperature_steps);
    write_1D_array_to_csv("U_block_hist_300_1000K.dat", energy_block_history, temperature_steps);

    write_1D_array_to_csv("C_hist_300_1000K.dat", capacity_history, temperature_steps);
    // write_1D_array_to_csv("C_hist_300_1000K.dat", capacity_corr_history, temperature_steps);
    // write_1D_array_to_csv("C_hist_300_1000K.dat", capacity_block_history, temperature_steps);

    write_1D_array_to_csv("long_hist_300_1000K.dat", long_order_history, temperature_steps);
    write_1D_array_to_csv("long_corr_hist_300_1000K.dat", long_order_corr_history, temperature_steps);
    write_1D_array_to_csv("long_block_hist_300_1000K.dat", long_order_history, temperature_steps);

    write_1D_array_to_csv("short_hist_300_1000K.dat", short_order_history, temperature_steps);
    write_1D_array_to_csv("short_corr_hist_300_1000K.dat", short_order_corr_history, temperature_steps);
    write_1D_array_to_csv("short_block_hist_300_1000K.dat", short_order_block_history, temperature_steps);

    // Free allocated memory
    free(temperature);
    free(energy_history);
    free(energy_corr_history);
    free(energy_block_history);

    free(capacity_history);
    free(capacity_corr_history);
    free(capacity_block_history);

    free(long_order_history);
    free(long_order_corr_history);
    free(long_order_block_history);

    free(short_order_history);
    free(short_order_corr_history);
    free(short_order_block_history);
}


double 
calculate_stat_ineff_autocorr(double *data, int data_length) {
    // Calculate the statistical inefficiency using the autocorrelation function.
    for (int k = 200; k < 75000; k++) {
        if (autocorrelation(data, data_length, k) < 0.125) {
            return (double) k;
        }
    }
    return 75000.;
}

double 
calculate_stat_ineff_block(double *data, int data_length) {
    // Calculate the statistical inefficiency using the block average method. 
    int block_length_min = 90000;
    int block_length_max = 100000;
    double statistical_inefficiency_block[block_length_max-block_length_min];
    for (int i = 0; i < block_length_max-block_length_min; i++) {
        statistical_inefficiency_block[i] = block_average(data, data_length, (i+block_length_min));
    }
    return average(statistical_inefficiency_block, block_length_max-block_length_min);
}
