#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>

#include "tools.h"
#include "H3_utils.h"

#define PI 3.141592653589

void task1();
void task2();

int
run(
    int argc,
    char *argv[]
   )
{
    task2();
    return 0;
}

void 
task1()
{
    // Simulation parameters
    int numb_walker_start = 200;
    double numb_walker_current = 200.;
    int time_steps = (int) 2e5; 
    int burnin_time_steps = (int) 5e3;
    double delta_t = 0.02;
    double average_E_T = 0.5;
    double E_T = 0.5;
    double damping_parameter = 0.5;
    
    // Array to save simulation results
    int total_pos_samples = 0;
    double *walker_positions = malloc(sizeof(double) * numb_walker_start * time_steps * 2);
    double *E_T_history = malloc(sizeof(double)*time_steps);
    double *numb_walker_history = malloc(sizeof(double)*time_steps);

    // Initializing GSL RNG
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    int seed = 42;
    gsl_rng_set(r, seed); 

    // Initialize walkers
    double *walker_list = malloc(numb_walker_start*15*sizeof(double));  // Creating larger array than necessary to be sure that an increase in # of walkers doesn't cause memory failure.
    double *equidistant_positions = linspace(-5., 5., numb_walker_start);
    for(int i = 0; i < numb_walker_start; i++) {
        walker_list[i] = equidistant_positions[i];
    }
    free(equidistant_positions);

    // Monte Carlo Sampling
    for(int i = 0; i < time_steps; i++) {
        // Print progress
        if(i % (time_steps/100) == 0){
            char chars[] = {'-', '\\', '|', '/'};
            printf(" %c | Progress: %d%% | i: %d \r", chars[(i % 99) % 4], (i * 100) / time_steps, i);
            fflush(stdout);
        }

        // Spreading of the wavefunction (walkers) due to kinetic energy.
        for(int j = 0; j < numb_walker_current; j++) {
            walker_list[j] += gsl_ran_gaussian(r, 1.0) * pow(delta_t, 1./2.);
        }

        // Life/Death process of walkers with weight function
        for(int j = 0; j < numb_walker_current; j++) {
            int m = numberOfWalkers(walker_list[j], delta_t, E_T, pot_task1, r);

            if(m == 0 && numb_walker_current > 1.5) {
                remove_element_at_index(walker_list, j, numb_walker_current);
                numb_walker_current -= 1;
                --j;
            } else if (m > 1 && numb_walker_current < numb_walker_start*15 - 2) {
                for(int k = 0; k < (m-1); k++){
                    walker_list[(int)numb_walker_current + k] = walker_list[j];
                }
                numb_walker_current += m - 1;
            }
        }

        // Update average E_T
        E_T = average_E_T - damping_parameter * log(numb_walker_current/numb_walker_start);
        if(i < burnin_time_steps) {
            // Updating average E_T with the burnin steps.
            average_E_T = calculate_average_E_T(i + 1, E_T, average_E_T);
        } else {
            // Updating average E_T without the burnin steps.
            average_E_T = calculate_average_E_T(i - burnin_time_steps + 1, E_T, average_E_T);
        }
        
        // Save the obtained values for E_T and number of walkers
        E_T_history[i] = average_E_T;
        numb_walker_history[i] = numb_walker_current;

        // Save the walkers positions
        for(int j = 0; j < numb_walker_current; j++) {
            walker_positions[total_pos_samples + j] = walker_list[j];
        }
        total_pos_samples += numb_walker_current;
    }

    //Saving objects to .dat files
    write_1D_array_to_csv("./data/task1/E_T_task1.dat", E_T_history, time_steps);
    write_1D_array_to_csv("./data/task1/NumbWalkers_task1.dat", numb_walker_history, time_steps);
    write_1D_array_to_csv("./data/task1/walker_pos_task1.dat", walker_positions, total_pos_samples);

    // free objects from memory
    gsl_rng_free(r);
    free(walker_list);
    free(E_T_history);
    free(numb_walker_history);
    free(walker_positions);
}

void 
task2()
{
    // Spinning process tracker
    int spin_tracker = 0;

    // Simulation parameters
    int numb_walker_start = 1000;
    double numb_walker_current = 1000.;
    int time_steps = (int) 5e5; 
    int burnin_time_steps = (int) 1e4;
    double delta_t = 0.01;
    double average_E_T = -3;
    double E_T = -3;
    double damping_parameter = 0.5;
    
    // Array to save simulation results
    double *E_T_history = malloc(sizeof(double)*time_steps);
    double *numb_walker_history = malloc(sizeof(double)*time_steps);

    // Initializing GSL RNG
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    int seed = 45;
    gsl_rng_set(r, seed); 

    // Initialize walkers
    double **walker_list = create_2D_array(numb_walker_start*500, 6);  // Creating larger array than necessary to be sure that an increase in # of walkers doesn't cause memory failure.
    for(int i = 0; i < numb_walker_start; i++) {
        double r1 = 0.7 + gsl_rng_uniform(r);
        double theta1 = acos(2 * gsl_rng_uniform(r) - 1);
        double phi1 = 2 * PI * gsl_rng_uniform(r);

        double r2 = 0.7 + gsl_rng_uniform(r);
        double theta2 = acos(2 * gsl_rng_uniform(r) - 1);
        double phi2 = 2 * PI * gsl_rng_uniform(r);

        walker_list[i][0] = r1 * sin(theta1) * cos(phi1);
        walker_list[i][1] = r1 * sin(theta1) * sin(phi1);
        walker_list[i][2] = r1 * cos(theta1);

        walker_list[i][3] = r2 * sin(theta2) * cos(phi2);
        walker_list[i][4] = r2 * sin(theta2) * sin(phi2);
        walker_list[i][5] = r2 * cos(theta2);
    }

    // Monte Carlo Sampling
    for(int i = 0; i < time_steps; i++) {
        // Print progress
        if(i % (time_steps/100) == 0){
            char chars[] = {'-', '\\', '|', '/'};
            printf("  %c | Progress: %d%% | i: %d \n", chars[spin_tracker % 4], (i * 100) / time_steps, i);

            spin_tracker++;
            fflush(stdout);
        }

        // Spreading of the wavefunction (walkers) due to kinetic energy.
        for(int j = 0; j < numb_walker_current; j++) {
            for(int k = 0; k < 6; k++) {
                walker_list[j][k] += gsl_ran_gaussian(r, 1.0) * pow(delta_t, 1./2.);
            }
        }

        // Life/Death process of walkers with weight function
        int *mult = malloc((int)numb_walker_current * sizeof(int));
        int total_new = 0;
        for(int j = 0; j < numb_walker_current; j++){
            mult[j] = numberOfWalkers2(walker_list[j], delta_t, E_T, pot_task2, r);
            total_new += mult[j];
        }

        double **tmp = create_2D_array(total_new, 6);
        int idx = 0;
        for(int j = 0; j < numb_walker_current; j++){
            for(int copy = 0; copy < mult[j]; copy++){
                for(int k=0; k<6; k++){
                    tmp[idx][k] = walker_list[j][k];
                }
                idx++;
            }
        }
        free(mult);

        if(total_new > 3000){
            for(int i=0; i < total_new; i++){  // Fisher–Yates shuffle
                int s = i + (int)(gsl_rng_uniform(r)*(total_new - i));
                double *swap = tmp[i]; tmp[i] = tmp[s]; tmp[s] = swap;
            }
            total_new = 3000;
        }
        for(int i=0; i<total_new; i++){
            for(int k=0; k<6; k++){
                walker_list[i][k] = tmp[i][k];
            }
        }
        destroy_2D_array(tmp, idx); // idx was the actual array size
        numb_walker_current = total_new;

        // Update average E_T
        E_T = average_E_T - damping_parameter * log(numb_walker_current/numb_walker_start);
        if(i < burnin_time_steps) {
            // Updating average E_T with the burnin steps.
            average_E_T = calculate_average_E_T(i + 1, E_T, average_E_T);
        } else {
            // Updating average E_T without the burnin steps.
            average_E_T = calculate_average_E_T(i - burnin_time_steps + 1, E_T, average_E_T);
        }
        
        // Save the obtained values for E_T and number of walkers
        E_T_history[i] = average_E_T;
        numb_walker_history[i] = numb_walker_current;
    }

    //Saving objects to .dat files
    write_1D_array_to_csv("./data/task2/E_T_task2.dat", E_T_history, time_steps);
    write_1D_array_to_csv("./data/task2/NumbWalkers_task2.dat", numb_walker_history, time_steps);

    // free objects from memory
    gsl_rng_free(r);
    destroy_2D_array(walker_list, numb_walker_start*500);
    free(E_T_history);
    free(numb_walker_history);
    //free(walker_positions);
}


