#pragma once
#include <stdio.h>

/* **********************************************
 *
 * Add uniform noise to a 2D array
 * The array elements after transformation is: 
 * The size of the noise is: +-noise. 
 *
 * **********************************************/
void
add_noise_2D_array(
    double **array, 
    double noise, 
    int rows, 
    int cols);


/* **********************************************
 *
 * Solves the equations of motions 
 * using potential.c and updates the 
 * positions using the verlat algorithm.
 *
 * **********************************************/
void velocity_verlet_one_step(
    double **positions, 
    double **velocities,
    double **forces, 
    double *virals,
    double *potential,
    double mass,
    double timestep,
    double cell_length,
    int nbr_atoms);


/* **********************************************
 *
 * Calculates the kintetic energy of n atoms
 * moving in 3D.
 * 
 * velocities on the form n_atoms x 3.
 *
 * **********************************************/
double 
calculate_kinetic_energy(
    double **velocities, 
    double mass, 
    int n_atoms);