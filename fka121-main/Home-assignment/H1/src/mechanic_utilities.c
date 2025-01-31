#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "mechanic_utilities.h"
#include "potential.h"

void
add_noise_2D_array(
    double **array, 
    double noise, 
    int rows, 
    int cols)     
{
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    int seed = 41;
    gsl_rng_set(r, seed); 

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            array[i][j] += (gsl_rng_uniform(r) - 0.5) * 2 * noise;
        }
    }


    gsl_rng_free (r);
}

void 
velocity_verlet_one_step(
    double **positions, 
    double **velocities,
    double **forces, 
    double *virals,
    double *potential,
    double mass,
    double timestep,
    double cell_length,
    int nbr_atoms)
{
    for (int i = 0; i < nbr_atoms; i++) {
        for(int j = 0; j < 3; j++) {
            positions[i][j] += velocities[i][j] * timestep + pow(timestep,2) * forces[i][j]/(mass * 2);
            velocities[i][j] +=  forces[i][j]/mass * timestep/2;
        }
    }
    
    calculate(potential, virals, forces, positions, cell_length, nbr_atoms);

    for (int i = 0; i < nbr_atoms; i++) {
        for(int j = 0; j < 3; j++) {
            velocities[i][j] +=  forces[i][j]/mass * timestep/2;
        }
    }
}

double 
calculate_kinetic_energy(
    double **velocities, 
    double mass, 
    int n_atoms)
{
    double kinetic_energy = 0.;
    for (int i = 0; i < n_atoms; i++) {
        for (int j = 0; j < 3; j++) {
            kinetic_energy += pow(velocities[i][j],2)*mass/2;            
        }
    } 
    return kinetic_energy;
}