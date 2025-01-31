#pragma once
#include <stdio.h>

#include <gsl/gsl_rng.h>

/* **********************************************
 *
 * Boltzmann constant, used in energy calculations. eV / K
 *
 * **********************************************/
extern const double k_B;

/* **********************************************
 *
 * Array of offsets representing the 8 corners of 
 * a cube relative to its center in a BCC lattice.
 *
 * **********************************************/
extern const int corner_offsets[8][3];

/* **********************************************
 *
 * Creates a 3D lattice initialized with a given atom type.
 *
 * lattice_size: Dimension of the lattice (cubic).
 * atom_type: The type of atom to initialize the lattice with.
 * 
 * Returns: Pointer to the 3D lattice array.
 *
 * **********************************************/
double ***create_cold_lattice(
    int lattice_size, 
    double atom_type);

/* **********************************************
 *
 * Calculates the total energy of a BCC lattice system.
 *
 * lattice_A: Pointer to the first lattice array.
 * lattice_B: Pointer to the second lattice array.
 * lattice_size: Dimension of the lattice (cubic).
 * E_AA: Interaction energy between two A atoms.
 * E_BB: Interaction energy between two B atoms.
 * E_AB: Interaction energy between an A and B atom.
 * 
 * Returns: Total energy of the system.
 *
 * **********************************************/
double calculate_BCC_lattice_energy(
    double ***lattice_A,
    double ***lattice_B,
    int lattice_size,
    double E_AA,
    double E_BB,
    double E_AB
);

/* **********************************************
 *
 * Calculates the long-range order parameter for a lattice.
 *
 * lattice_A: Pointer to the lattice array.
 * lattice_size: Dimension of the lattice (cubic).
 * 
 * Returns: Long-range order parameter.
 *
 * **********************************************/
double calculate_long_range_order(
    double ***lattice_A,
    int lattice_size
);

/* **********************************************
 *
 * Calculates the short-range order parameter for a lattice system.
 *
 * lattice_A: Pointer to the first lattice array.
 * lattice_B: Pointer to the second lattice array.
 * lattice_size: Dimension of the lattice (cubic).
 * 
 * Returns: Short-range order parameter.
 *
 * **********************************************/
double calculate_short_range_order(
    double ***lattice_A, 
    double ***lattice_B,
    int lattice_size
);

/* **********************************************
 *
 * Performs one step of the Metropolis Monte Carlo algorithm.
 *
 * lattice_A: Pointer to the first lattice array.
 * lattice_B: Pointer to the second lattice array.
 * lattice_size: Dimension of the lattice (cubic).
 * E_AA: Interaction energy between two A atoms.
 * E_BB: Interaction energy between two B atoms.
 * E_AB: Interaction energy between an A and B atom.
 * temp: Temperature of the system.
 * energy: Current energy of the system.
 * r: Pointer to the GSL random number generator.
 * 
 * Returns: 0 if no swap was performed, 1 if a swap was accepted,
 *          2 if a swap was rejected.
 *
 * **********************************************/
int metropolis_one_step(
    double ***lattice_A,
    double ***lattice_B,
    int lattice_size,
    double E_AA,
    double E_BB,
    double E_AB, 
    double temp,
    double *energy,
    gsl_rng *r
);

/* **********************************************
 *
 * Performs one step of the Metropolis Monte Carlo algorithm.
 *
 * Optimized by chatGPT.
 * 
 * lattice_A: Pointer to the first lattice array.
 * lattice_B: Pointer to the second lattice array.
 * lattice_size: Dimension of the lattice (cubic).
 * E_AA: Interaction energy between two A atoms.
 * E_BB: Interaction energy between two B atoms.
 * E_AB: Interaction energy between an A and B atom.
 * temp: Temperature of the system.
 * energy: Current energy of the system.
 * r: Pointer to the GSL random number generator.
 * 
 * Returns: 0 if no swap was performed, 1 if a swap was accepted,
 *          2 if a swap was rejected.
 *
 * **********************************************/
int 
metropolis_one_step_optimized(
    double ***lattice_A,
    double ***lattice_B,
    int lattice_size,
    double E_AA,
    double E_BB,
    double E_AB, 
    double temperature,
    double *energy,
    gsl_rng *r);