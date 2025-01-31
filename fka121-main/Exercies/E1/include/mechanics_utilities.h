#pragma once
#include <stdio.h>

// Function to calculate acceleration
// Parameters:
//  - accelerations: Array to store the calculated accelerations
//  - positions: Array of positions
//  - masses: Array of masses
//  - kappa: Spring constant
void calculate_acceleration(double *accelerations, double *positions, double *masses, double kappa);

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
void calculate_acceleration_nonlinear(double *accelerations, double *positions, double alpha, const unsigned int N);

// Function to calculate potential energy
// Parameters:
//  - positions: Array of positions
//  - kappa: Spring constant
// Returns:
//  - The calculated potential energy
double calculate_potential_energy(double *positions, double kappa);

// Function to calculate kinetic energy
// Parameters:
//  - velocities: Array of velocities
//  - masses: Array of masses
// Returns:
//  - The calculated kinetic energy
double calculate_kinetic_energy(double *velocities, double *masses);

// Function to perform one step of the velocity Verlet algorithm
// Parameters:
//  - accelerations: Array of accelerations
//  - positions: Array of positions
//  - velocities: Array of velocities
//  - masses: Array of masses
//  - kappa: Spring constant
//  - timestep: The time step for the integration
void velocity_verlet_one_step(double *accelerations, double *positions, double *velocities,
                              double *masses, double kappa, double timestep);

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
                              double alpha, double timestep, const unsigned int N);

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
void transform_to_normal_modes(double *positions, double *Q, const unsigned int N);

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
void calculate_normal_mode_energies(double *energies, double *positions, double *velocities, const unsigned int N);
