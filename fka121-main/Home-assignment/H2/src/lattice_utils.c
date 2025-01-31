#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "tools.h"

const double k_B = 8.617333262e-5;

const int corner_offsets[8][3] = {
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
        {1, 1, 0},
        {1, 0, 1},
        {0, 1, 1},
        {1, 1, 1}
};

double 
***create_cold_lattice(int lattice_size, double atom_type) {
    double ***lattice = create_3D_array(lattice_size, lattice_size, lattice_size);

    for(int i = 0; i < lattice_size; i++) {
        for (int j = 0; j < lattice_size; j++) {
            for (int k = 0; k < lattice_size; k++) {
                lattice[i][j][k] = atom_type;
            }
        }
    }

    return lattice;
}


double
calculate_BCC_lattice_energy(
    double ***lattice_A,
    double ***lattice_B,
    int lattice_size,
    double E_AA,
    double E_BB,
    double E_AB)
{
    double total_energy = 0.0;

    /* The BCC lattice can be viewed as two interpenetrating simple cubic lattices: 
       Lattice_B nodes occupy positions shifted by half a cell relative to Lattice_A nodes.
       Each unit cell in the BCC lattice consists of a central atom (from Lattice_B) and
       eight surrounding corner atoms (from Lattice_A).
       
       The following offsets define the relative positions of these eight corner atoms 
       with respect to the central atom at lattice_B[i][j][k]. 
       These offsets represent all 8 corner points of a cube around the central point.
       
       If the central atom is at (i, j, k) in Lattice_B, then the corner atoms are found at:
       (i + 0, j + 0, k + 0),
       (i + 1, j + 0, k + 0),
       (i + 0, j + 1, k + 0),
       (i + 0, j + 0, k + 1),
       (i + 1, j + 1, k + 0),
       (i + 1, j + 0, k + 1),
       (i + 0, j + 1, k + 1),
       (i + 1, j + 1, k + 1).
       
       These 8 points correspond exactly to the corners of a cube whose center is 
       occupied by the atom from Lattice_B at (i, j, k). 
    */

    // Iterate over all positions in Lattice_B (the "center" atoms of the BCC unit cells).
    for (int i = 0; i < lattice_size; i++) {
        for (int j = 0; j < lattice_size; j++) {
            for (int k = 0; k < lattice_size; k++) {
                // Get the type of the central atom (1 or 2)
                int center_type = (int) lattice_B[i][j][k];

                // Loop over all eight corners of the cube around the central atom
                for (int n = 0; n < 8; n++) {
                    int ia = i + corner_offsets[n][0];
                    int ja = j + corner_offsets[n][1];
                    int ka = k + corner_offsets[n][2];

                    // Check if the corner index is within the lattice bounds and else apply periodic boundry conditions.
                    if (ia >= lattice_size) {
                        ia = 0;
                    } 
                    if (ja >= lattice_size) {
                        ja = 0;
                    } 
                    if (ka >= lattice_size) {
                        ka = 0;
                    }

                    // Get the type of the corner atom
                    int corner_type = (int) lattice_A[ia][ja][ka];

                    // Add the corresponding interaction energy:
                    // E_AA if both are type 1 (A)
                    // E_BB if both are type 2 (B)
                    // E_AB otherwise
                    if (center_type == 1 && corner_type == 1) total_energy += E_AA;
                    else if (center_type == 2 && corner_type == 2) total_energy += E_BB;
                    else total_energy += E_AB;

                    // printf("total energy %lf \n", total_energy);
                }
            }
        }
    }
    return total_energy;
}

double
calculate_long_range_order(
    double ***lattice_A,
    int lattice_size)
{
    double N_A = 0;
    for(int i = 0; i < lattice_size; i++) {
        for (int j = 0; j < lattice_size; j++) {
            for (int k = 0; k < lattice_size; k++) {
                N_A += lattice_A[i][j][k] - 1;
            }
        }
    }
    return N_A * 2 / pow((double) lattice_size, 3) - 1;
}

double
calculate_short_range_order(
    double ***lattice_A, 
    double ***lattice_B,
    int lattice_size) 
{
    double q = 0;
    for (int i = 0; i < lattice_size; i++) {
        for (int j = 0; j < lattice_size; j++) {
            for (int k = 0; k < lattice_size; k++) {
                // Get the type of the central atom (1 or 2)
                double center_type = lattice_B[i][j][k];

                // Loop over all eight corners of the cube around the central atom
                for (int n = 0; n < 8; n++) {
                    int ia = i + corner_offsets[n][0];
                    int ja = j + corner_offsets[n][1];
                    int ka = k + corner_offsets[n][2];

                    // Check if the corner index is within the lattice bounds and else apply periodic boundry conditions.
                    if (ia >= lattice_size) {
                        ia = 0;
                    } 
                    if (ja >= lattice_size) {
                        ja = 0;
                    } 
                    if (ka >= lattice_size) {
                        ka = 0;
                    }

                    // Get the type of the corner atom
                    double corner_type = lattice_A[ia][ja][ka];

                    // Increment q with 1 if the neighbours are of different type
                    if ((center_type == 1. && corner_type == 2.) || (center_type == 2. && corner_type == 1.)) {
                        q += 1.;
                    }
                }
            }
        }
    }

    double N = pow((double) lattice_size, 3);
    return 1 / (4 * N) * (q - 4*N); 
}


int 
metropolis_one_step(
    double ***lattice_A,
    double ***lattice_B,
    int lattice_size,
    double E_AA,
    double E_BB,
    double E_AB, 
    double temperature,
    double *energy,
    gsl_rng *r) 
{
    int ib = gsl_rng_uniform_int(r, lattice_size);
    int jb = gsl_rng_uniform_int(r, lattice_size);
    int kb = gsl_rng_uniform_int(r, lattice_size);

    int ia = gsl_rng_uniform_int(r, lattice_size);
    int ja = gsl_rng_uniform_int(r, lattice_size);
    int ka = gsl_rng_uniform_int(r, lattice_size);

    if ((int) lattice_A[ia][ja][ka] == (int) lattice_B[ib][jb][kb]) {
        return 0;
    } else {
        double temp = lattice_A[ia][ja][ka];
        lattice_A[ia][ja][ka] = lattice_B[ib][jb][kb];
        lattice_B[ib][jb][kb] = temp;

        double new_energy = calculate_BCC_lattice_energy(lattice_A, lattice_B, lattice_size, E_AA, E_BB, E_AB); 

        double step_probability = exp(-(new_energy - *energy)/(k_B * temperature));
        double step_draw = gsl_rng_uniform(r);

        if (step_probability >= step_draw) {
            *energy = new_energy;
            return 1;
        } else {
            double temp = lattice_A[ia][ja][ka];
            lattice_A[ia][ja][ka] = lattice_B[ib][jb][kb];
            lattice_B[ib][jb][kb] = temp;
            return 2;
        }
    }
}

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
    gsl_rng *r) 
{
    // Select random positions in lattice A and B
    int ia = gsl_rng_uniform_int(r, lattice_size);
    int ja = gsl_rng_uniform_int(r, lattice_size);
    int ka = gsl_rng_uniform_int(r, lattice_size);

    int ib = gsl_rng_uniform_int(r, lattice_size);
    int jb = gsl_rng_uniform_int(r, lattice_size);
    int kb = gsl_rng_uniform_int(r, lattice_size);

    // Check if the selected atoms are of the same type
    if ((int) lattice_A[ia][ja][ka] == (int) lattice_B[ib][jb][kb]) {
        return 0; // No swap performed
    } else {
        // Compute energy before swap for both atoms
        // For simplicity, assume interactions are local and only depend on neighbors
        // This requires implementation of a local energy calculation function
        double delta_E = 0.0;

        // Swap the atoms
        double temp = lattice_A[ia][ja][ka];
        lattice_A[ia][ja][ka] = lattice_B[ib][jb][kb];
        lattice_B[ib][jb][kb] = temp;

        // Compute energy after swap
        double new_energy = calculate_BCC_lattice_energy(lattice_A, lattice_B, lattice_size, E_AA, E_BB, E_AB); 

        delta_E = new_energy - (*energy);

        // Metropolis criterion
        if (delta_E <= 0 || gsl_rng_uniform(r) < exp(-delta_E / (k_B * temperature))) {
            *energy = new_energy;
            return 1; // Swap accepted
        } else {
            // Revert the swap
            temp = lattice_A[ia][ja][ka];
            lattice_A[ia][ja][ka] = lattice_B[ib][jb][kb];
            lattice_B[ib][jb][kb] = temp;
            return 2; // Swap rejected
        }
    }
}



    