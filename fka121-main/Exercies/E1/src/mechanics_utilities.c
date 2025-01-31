#include <math.h>

#include "/home/felixpersson/kompfys/fka121-main/Exercies/include/mechanics_utilities.h"
#include <gsl/gsl_math.h>


/*
 * Calculate the acceleration.
 *
 * All vectors are assumed to be of length 3
 *
 * Parameters
 * ----------
 *  accelerations - Vector where accelerations are to be written
 *  positions - Vector of positions
 *  masses - Vector of masses
 *  kappa - Spring constant
 *
*/
void calculate_acceleration(double *accelerations, double *positions, double *masses, double kappa)
{
    const unsigned int N = 3;  // There are three particles
    for (int i = 0; i < N; i++) {
        if (i == 0) {
            accelerations[i] = kappa/masses[i] * (positions[i + 1] - positions[i]);
        } else if (i == (N - 1)) {
            accelerations[i] = kappa/masses[i] * (positions[i - 1] - positions[i]);
        } else {
            accelerations[i] = kappa/masses[i] * (positions[i - 1] - 2 * positions[i] + positions[i + 1]);
        }
    }
}

/**
 * Calculate the acceleration with a non-linear term.
 *
 * Assume that m = kappa = 1
 *
 * Parameters
 * ----------
 *  accelerations - Vector where accelerations are to be written
 *  positions - Vector of positions
 *  alpha - Anharmonicity constant
 *  N - Number of atoms
 *
*/
void calculate_acceleration_nonlinear(double *accelerations, double *positions, double alpha, const unsigned int N)
{
    for (int i = 0; i < N; i++) {
        if (i == 0) {
            accelerations[i] = (positions[i + 1] - 2*positions[i]) + alpha * (pow((positions[i + 1] - positions[i]),2) - pow(positions[i],2));
        } else if (i == (N - 1)) {
            accelerations[i] = (positions[i - 1] - 2*positions[i]) + alpha * (pow((positions[i]),2) - pow(positions[i] - positions[i-1],2));;
        } else {
            accelerations[i] = (positions[i - 1] - 2 * positions[i] + positions[i + 1]) + alpha * (pow((positions[i + 1] - positions[i]),2) - pow(positions[i] - positions[i - 1],2));
        }
    }
}

/**
 * Calculate the potential energy
 *
 * All vectors are assumed to be of length 3
 *
 * Parameters
 * ----------
 *  positions - Vector of positions
 *  kappa - Spring constant
 * Returns
 * -------
 *  Potential energy
 *
*/
double calculate_potential_energy(double *positions, double kappa)
{
    const unsigned int N = 3;  // There are three particles
    double potential = 0.;
    for (int i = 0; i<(N-1); i++) {
        potential += kappa/2 * pow((positions[i] - positions[i + 1]), 2);
    }
    return potential;
}

/**
 * Calculate the kinetic energy
 *
 * All vectors are assumed to be of length 3
 *
 * Parameters
 * ----------
 *  velocities - Vector of velocities
 *  masses - Vector of masses
 *
 * Returns
 * -------
 *  Kinetic energy
 *
*/
double calculate_kinetic_energy(double *velocities, double *masses)
{
    const unsigned int N = 3;  // There are three particles
    double kinetic_energy = 0.;
    for (int i = 0; i < N; i++) {
        kinetic_energy += pow(velocities[i],2)*masses[i]/2;
    } 
    return kinetic_energy;
}

/**
 * Perform one velocity Verlet step
 *
 * All vectors are assumed to be of length 3
 *
 * Parameters
 * ----------
 *  accelerations - Vector of accelerations
 *  positions - Vector of positions
 *  velocities - Vector of velocities
 *  masses - Vector of masses
 *  kappa - Spring constant
 *  timestep - Time step
 *
*/
void velocity_verlet_one_step(double *accelerations, double *positions, double *velocities,
                              double *masses, double kappa, double timestep)
{
    const unsigned int N = 3;  // There are three particles
    
    for (int i = 0; i < N; i++) {
        positions[i] += velocities[i] * timestep + pow(timestep,2) * accelerations[i]/2;
        velocities[i] +=  accelerations[i] * timestep/2;
    }

    for (int i = 0; i < N; i++) {
        calculate_acceleration(accelerations, positions, masses, kappa);
        velocities[i] = velocities[i] + timestep * accelerations[i]/2;
    }
}

/**
 * Perform one velocity Verlet step with nonlinear force
 *
 * Assume m = kappa = 1
 * 
 * Parameters
 * ----------
 *  accelerations - Vector of accelerations
 *  positions - Vector of positions
 *  velocities - Vector of velocities
 *  alpha - Anharmonicity constant
 *  timestep - Time step
 *  N - Number of atoms
 *
*/
void velocity_verlet_one_step_nonlinear(double *accelerations, double *positions, double *velocities,
                              double alpha, double timestep, const unsigned int N)
{
    for (int i = 0; i < N; i++) {
        positions[i] += velocities[i] * timestep + pow(timestep,2) * accelerations[i]/2;
        velocities[i] +=  accelerations[i] * timestep/2;
    }

    for (int i = 0; i < N; i++) {
        calculate_acceleration_nonlinear(accelerations, positions, alpha, N);
        velocities[i] = velocities[i] + timestep * accelerations[i]/2;
    }
}

/**
 * Transform to normal modes.
 *
 * Parameters
 * ----------
 *  positions - Vector of positions
 *  Q - Vector of Q-coordinates
 *  N - Number of particles (length of vectors)
 *
*/
void transform_to_normal_modes(double *positions, double *Q, const unsigned int N)
{
    for (int k = 0; k < N; k++) {
        Q[k] = 0;
        for (int i = 0; i < N; i++) {
            Q[k] += (double) sqrt(2.0 / (N + 1.0)) * positions[i] * sin((i + 1)  * (k+1) * M_PI/(N+1.0));
        }
    }
}


/**
 * Calculate the normal mode energies
 *
 * Parameters
 * ----------
 *  energies - Vector where energies will be written
 *  positions - Vector of positions
 *  velocities - Vector of positions
 *  N - Number of particles (length of vectors)
 *
*/
void calculate_normal_mode_energies(double *energies, double *positions, double *velocities, const unsigned int N)
{
    double P[N];
    double Q[N];

    transform_to_normal_modes(velocities, P, N);
    transform_to_normal_modes(positions, Q, N);

    double omega[N];
    for (int k = 0; k < N; k++) {
        omega[k] = 2 * sin((k + 1.0) * M_PI / (2 * (N + 1.0))); 
    }

    for (int k = 0; k < N; k++) {
        energies[k] = (pow(P[k], 2) + pow(omega[k] * Q[k], 2))/2.0;
    }
}

