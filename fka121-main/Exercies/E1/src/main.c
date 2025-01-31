#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "/home/felixpersson/kompfys/fka121-main/Exercies/include/tools.h"
#include "/home/felixpersson/kompfys/fka121-main/Exercies/include/mechanics_utilities.h"

void task_1();
void simulate_Verlat_CO2(
    double *positions, 
    double *velocitys,
    double *accelerations,
    double *masses, 
    double kappa, 
    double time_step, 
    int N, 
    double **positions_history,
    double *kinetic_energy,
    double *potential_energy);

void task_3();
void simulate_Verlat_harmonic(
    double *eigen_energies,
    double *positions, 
    double *velocitys,
    double *accelerations,
    double time_step,
    double alpha, 
    int iterations,
    int N,  
    double **energy_history);

void task_4();

void simulate_long_Verlat_harmonic(
    double *eigen_energies,
    double *positions, 
    double *velocitys,
    double *accelerations,
    double time_step,
    double alpha, 
    int iterations,
    int N,
    int M,  
    double **energy_history);

int
main()
{
    task_4();

    return 0;
}


void task_1() {
    double positions[3] = {0.01, 0.005, -0.005};
    double velocitys[3] = {0,0,0};
    double accelerations[3] = {0,0,0};
    double masses[3] = {0.001658, 0.00124479, 0.001658};
    double kappa = 99.86;

    int iterations = 500;
    double time_step = 0.005;

    double **positions_history = create_2D_array(iterations, 3);
    double *kinetic_energy = calloc(iterations, sizeof(double)); 
    double *potential_energy = calloc(iterations, sizeof(double)); 

    simulate_Verlat_CO2(
        positions, 
        velocitys, 
        accelerations, 
        masses, 
        kappa, 
        time_step, 
        iterations, 
        positions_history,
        kinetic_energy, 
        potential_energy);
}

void simulate_Verlat_CO2(
    double *positions, 
    double *velocitys,
    double *accelerations,
    double *masses, 
    double kappa, 
    double time_step, 
    int N, 
    double **positions_history,
    double *kinetic_energy,
    double *potential_energy) 
{
    for (int i = 0; i < 3; i++) {
        positions_history[0][i] = positions[i];
    }

    for (int i = 0; i < N; i++) {
        velocity_verlet_one_step(accelerations, positions, velocitys, masses, kappa, time_step);
        for (int j = 0; j < 3; j++) {
            positions_history[i][j] = positions[j];
        }
        kinetic_energy[i] = calculate_kinetic_energy(velocitys, masses);
        potential_energy[i] = calculate_potential_energy(positions, kappa);
    }
    /*
    write_2D_array_to_csv("Positions.dat", positions_history, N, 3);
    write_1D_array_to_csv("Potential_Energy.dat", potential_energy, N);
    write_1D_array_to_csv("Kinetic_Energy.dat", kinetic_energy, N);
    */
    free(kinetic_energy);
    free(potential_energy);
    for (int i = 0; i < N; i++) {
        free(positions_history[i]);
    }
    free(positions_history);
}

void task_3() {
    int N = 32;
    double alpha = 0;

    double Q[32] = {0};
    double P[32] = {0};
    P[0] = sqrt(2. * N);

    double positions[32];
    double velocities[32];
    double accelerations[32] = {0};

    transform_to_normal_modes(Q, positions, N);
    transform_to_normal_modes(P, velocities, N);

    double eigen_energies[32];
    calculate_normal_mode_energies(eigen_energies, positions, velocities, N);

    int iterations = 250000;
    double time_step = 0.1;

    double **energy_history = create_2D_array(iterations, 5);

    simulate_Verlat_harmonic(
        eigen_energies,
        positions, 
        velocities,
        accelerations,
        time_step,
        alpha, 
        iterations,
        N,  
        energy_history);

    write_2D_array_to_csv("Normal_mode_Energies.dat", energy_history, iterations, 5);

    for (int i = 0; i < iterations; i++) {
        free(energy_history[i]);
    }
    free(energy_history);
}

void simulate_Verlat_harmonic(
    double *eigen_energies,
    double *positions, 
    double *velocitys,
    double *accelerations,
    double time_step,
    double alpha, 
    int iterations,
    int N,  
    double **energy_history) 
{
    for (int i = 0; i < iterations; i++) {
        for (int j = 0; j < 5; j++) {
            energy_history[i][j] = eigen_energies[j];
        }

        velocity_verlet_one_step_nonlinear(accelerations, positions, velocitys, alpha, time_step, N);
        calculate_normal_mode_energies(eigen_energies, positions, velocitys, N);

        if ((i+1) % 25000 == 0) {
            printf("Progress: %d %% \n", i*100/iterations);
            for (int k = 0; k < 5; k++) {
                printf("energy_history[%d][%d] = %lf \n", i,k,energy_history[i][k]);
            }
        }
    }
}

void task_4() {
    int N = 32;
    int M = 1000;
    double alpha = 0.16;

    double Q[32] = {0};
    double P[32] = {0};
    P[0] = sqrt(2. * N);

    double positions[32];
    double velocities[32];
    double accelerations[32] = {0};

    transform_to_normal_modes(Q, positions, N);
    transform_to_normal_modes(P, velocities, N);

    double eigen_energies[32];
    calculate_normal_mode_energies(eigen_energies, positions, velocities, N);

    int iterations = 100000000; //10^7 iterations and t_max = 10^6.
    double time_step = 0.1;

    double **energy_history = create_2D_array(iterations/M, 32);

    simulate_long_Verlat_harmonic(
        eigen_energies,
        positions, 
        velocities,
        accelerations,
        time_step,
        alpha, 
        iterations,
        N,  
        M,
        energy_history);

    write_2D_array_to_csv("Long_Normal_mode_Energies_alpha01.dat", energy_history, iterations/M, 32);

    for (int i = 0; i < iterations/M; i++) {
        free(energy_history[i]);
    }
    free(energy_history);
}

void simulate_long_Verlat_harmonic(
    double *eigen_energies,
    double *positions, 
    double *velocities,
    double *accelerations,
    double time_step,
    double alpha, 
    int iterations,
    int N,  // Number of normalmodes/Atoms 
    int M,  // Steps between each saving of the normal mode energies
    double **energy_history) 
{
    for (int i = 0; i < iterations/M; i++) {
        for(int j = 0; j < M; j++) {
            for (int k = 0; k < 32; k++) {
                energy_history[i][k] = eigen_energies[k];
            }
        }
        velocity_verlet_one_step_nonlinear(accelerations, positions, velocities, alpha, time_step, N);
        calculate_normal_mode_energies(eigen_energies, positions, velocities, N);

        if ((i+1) % 1000 == 0) {
            printf("%d: Progress: %d %% \n",i, (int) i*100*M/iterations);
        }
    }
}