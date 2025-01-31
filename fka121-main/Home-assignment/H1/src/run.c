#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"
#include "/home/felixpersson/kompfys/fka121-main/Home-assignment/H1/include/mechanic_utilities.h"

const double k_B = 8.617333262e-5; // eV/K 

//Function declarations
void task_1();
void task_2();
void task_3();
void task_4();

void simulate_Verlat_task_2(
    double **positions, 
    double **velocities,
    double **forces, 
    double *virals,
    double *potential,
    double mass,
    double cell_length,
    int nbr_atoms,
    int iterations, 
    double time_step,
    double *kinetic_history,
    double *potential_history);

void simulate_Verlat_task_3(
    double **positions, 
    double **velocities,
    double **forces, 
    double *virials,
    double *potential,
    double mass,
    double isothermal_compression,
    double cell_length,
    double T_eq, 
    double P_eq,
    int nbr_atoms,
    int iterations, 
    double time_step,
    double time_constant,
    double *kinetic_history,
    double *virial_history,
    double *cell_length_histroy, 
    double **positions_history);

void simulate_Verlat_task_4(
    double **positions, 
    double **velocities,
    double **forces, 
    double *virials,
    double *potential,
    double mass,
    double isothermal_compression,
    double cell_length,
    double T_eq_high,
    double T_eq_low,
    double P_eq,
    int nbr_atoms,
    int iterations_high_temp,
    int iterations_low_temp, 
    double time_step,
    double time_constant,
    double *kinetic_history,
    double *virial_history,
    double *cell_length_history, 
    double **positions_history,
    double **velocities_history);     

//Function definitions
int
run(
    int argc,
    char *argv[]
   )
{       
    task_4();
    return 0;
}


// --------------------------------- Task 1 -------------------------------
void task_1() {
    int n_atoms = 256;
    int num_consts = 9;
    int N = 4;

    double *lattice_consts = linspace(4, 4.1, num_consts);
    double E_pot[num_consts];
    
    for (int i = 0; i < num_consts; i++) {
        double **positions = create_2D_array(n_atoms, 3);
        init_fcc(positions, N, lattice_consts[i]);

        E_pot[i] = get_energy_AL(positions, N * lattice_consts[i], n_atoms);
        destroy_2D_array(positions, n_atoms);
    }

    printf("Energies: ");
    for (int i = 0; i < num_consts; i++) {
        printf("%lf, ", E_pot[i]);
    }
    printf("\nLattice constants: ");
    for (int i = 0; i < num_consts; i++) {
        printf("%lf, ", lattice_consts[i]);
    }
    printf("\n");

    free(lattice_consts);
}

// --------------------------------- Task 2 -------------------------------

void task_2() {
    int n_atoms = 256;
    int N = 4; // size of supercell in each dimension.
    double mass = 2.7963e-3; // Al mass.
    double lattice_const = 4.;
    double cell_length = lattice_const * N;
    
    double **positions = create_2D_array(n_atoms, 3);
    init_fcc(positions, N, lattice_const);
    add_noise_2D_array(positions, 0.065 * lattice_const, n_atoms, 3);

    double **velocities = create_2D_array(n_atoms, 3);
    for (int i = 0; i < n_atoms; i++) {
        for (int j = 0; j < 3; j++) {
            velocities[i][j] = 0;
        }
    }
    double **forces = create_2D_array(n_atoms, 3);

    double virals = 0;
    double potential = 0;

    double time_step = 0.015;
    int iterations = 1000;

    double kinetic_history[iterations];  
    double potential_history[iterations];
    
    simulate_Verlat_task_2(
        positions, 
        velocities,
        forces, 
        &virals,
        &potential,
        mass,
        cell_length,
        n_atoms,
        iterations, 
        time_step,
        kinetic_history,
        potential_history);

    write_1D_array_to_csv("kinetic_energy_task2_0015.dat", kinetic_history, iterations);
    write_1D_array_to_csv("potential_energy_task2_0015.dat", potential_history, iterations);

    destroy_2D_array(positions, n_atoms); 
    destroy_2D_array(velocities, n_atoms);
    destroy_2D_array(forces, n_atoms);    
}

void simulate_Verlat_task_2(
    double **positions, 
    double **velocities,
    double **forces, 
    double *virals,
    double *potential,
    double mass,
    double cell_length,
    int nbr_atoms,
    int iterations, 
    double time_step,
    double *kinetic_history,
    double *potential_history) 
{
    for (int i = 0; i < iterations; i++) {
        potential_history[i] = *potential;
        kinetic_history[i] = calculate_kinetic_energy(velocities, mass, nbr_atoms);

        velocity_verlet_one_step(
            positions, 
            velocities,
            forces, 
            virals,
            potential,
            mass,
            time_step,
            cell_length,
            nbr_atoms);
        
        if ((i+1) % (iterations/20) == 0) {
            printf("Progress: %3.1lf %%\n", (double) (i+1)/iterations * 100);
        } 
    }
}

// --------------------------------- Task 3 -------------------------------

void task_3() {
    int n_atoms = 256;
    int N = 4;  // size of supercell in each dimension.
    double mass = 2.7963e-3; // Al mass.
    double lattice_const = 4.;  // Å
    double cell_length = lattice_const * N;

    double bulk_modulus = 70e9 * 6.24e-12; //  eV/Å^3, https://pdf.sciencedirectassets.com/271577/1-s2.0-S0038109800X1455X/1-s2.0-S0038109801005178/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEG4aCXVzLWVhc3QtMSJHMEUCIGeGoeFuBWjXDXfcOs9PKu4Cd9b3xlOaBR%2FrxrrM4Pj%2FAiEA2YNEKeGAywx9m0KMcgc1po0En8SL2Qjz97TIx0nZ3JMqswUIFxAFGgwwNTkwMDM1NDY4NjUiDDVnZR9Zb86Gc5qqtSqQBfrR3NlbxhXrFueJyW6qkO5hTZ7%2BcCTuqPOeSNmFRI26k6aOA4KQzjtZhZi6LnAve3e1a9yN01UYFL%2F4K6nyXL69%2Bsr7%2B5v2pQtxTSq7Q%2BNJXT75D2W64aepYHduL%2BAMDYcLBntpB4D3AxxP3lFDmmmICuOxGILnn2ZWvr0jJPyDkgOcskXRNwriTPiQ%2BbLjJwruGji5HGsDAs7hXLbwDa7EwYVF6ATkyYZMrkOQJ9J027sk8oyAWBX9lHzfXqUm9NsEJDi1kJYEhLrWcHJO42AB613EpiRfO1wSq9rzJF4SmFI%2FE8g3YaYZOal%2B3aH%2BjOhi%2B1q84Bu5fkE%2BWWCUPYnkG8s0aeaTFDFy4O%2FuZ9MHlOHgp3ZOnES8oilMqbdtDM9cNSOqjM4sVB%2Fw5x18qoXNNYWdtex6ix519UEjQUtINd%2FRqOK9I0gpVsId0CIpKvvM5II3RIScmjaJoAMTyCCelxbHWZWWzbR%2BBHroa02MFxekGOiCy6cunDZeNCxAbd%2BXZdi%2BsXnW9lA8G%2F%2FTpn2M6Uac6x6ZtC7rZSzS0Z2%2BcJiQRufh%2BhiSmpVjexUsLXuOpUZbA8XYi32G8%2F2FpK%2FUBepP9jqtAKBUI6YWYkTU8zMxyj3cdTToX8Vn7%2BUc4ut4i6JI9Eq%2Bw8Hh1xkUimu7d5toqbuIDRISspyMGMkF3ml68j1%2FLsNUquT9KV7wi40osXXU0d%2Bs8%2FkLbjeY2jkRsDVrKduTdH40%2FsJlyMBgbV0KowN1SNP2IChjHsUGP%2Bkgb8o7goAHWpzs8MLyVx5jsGLx%2FvxerTTyg6UnhDNJc9Iidvutd4M1WLYDBD1C442HOUVsjYoxzn2JfrnemJ2lBdGvFkfR61bA46xZTsAtMIv4kboGOrEBGADuQ1%2BQlHCvtyEjbIRWPw3xVW30goSR0y0sljN%2BydcVTlnIKaAxN6pcnyzJSyrb%2FHOSDzj9rLx3CY2zXIjIaCBOcc3G7u5gcAfiLS4Vi4ku%2FGqnOaKcazoBglmz4Hj4Nwo4PsdkWyuufrWf5klKVslPYoUIb%2FCjWjI2N6CP7OKmXs0QAC03BQUP5RLJ5HAsyJhqsvR8jldO8FyQQTkChF61zubt%2BSbgT9cddC7tMkN3&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20241125T140307Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYQQGRCNVY%2F20241125%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=ae7fba687d5efc6af8658cbfb805c168673ba6e5061e7d82d0c69fb5aa909772&hash=3ec358775603ee744e95194c31b3e7c635c5ff9545db2af16ba4e34ea442ca9a&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0038109801005178&tid=spdf-314dbf91-09ad-4a81-b820-80c10daef50a&sid=2b46fa2c46952945dc4862b74f1295abb677gxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=15055e015e505d005101&rr=8e822e2bac0febd0&cc=se
    double isothermal_compression = 1/bulk_modulus; // Å^3/eV

    double T_eq = 500. + 273.15;     // Kelvin
    double P_eq = 1e5 * 6.24e-12;    //eV/Å^3
    
    double **positions = create_2D_array(n_atoms, 3);
    init_fcc(positions, N, lattice_const);
    add_noise_2D_array(positions, 0.065 * lattice_const, n_atoms, 3);

    double **velocities = create_2D_array(n_atoms, 3);
    for (int i = 0; i < n_atoms; i++) {
        for (int j = 0; j < 3; j++) {
            velocities[i][j] = 0;
        }
    }
    double **forces = create_2D_array(n_atoms, 3);

    double virials = 0;
    double potential = 0;

    double time_step = 0.001;
    double time_const = 0.5;
    int iterations = 8500;

    double kinetic_history[iterations];  
    double virial_history[iterations];
    double cell_length_history[iterations];
    double **position_history = create_2D_array(iterations, 3);
    
     simulate_Verlat_task_3(
        positions, 
        velocities,
        forces, 
        &virials,
        &potential,
        mass,
        isothermal_compression,
        cell_length,
        T_eq,
        P_eq,
        n_atoms,
        iterations, 
        time_step,
        time_const,
        kinetic_history,
        virial_history,
        cell_length_history, 
        position_history);

    write_1D_array_to_csv("virial_history_T500_P1.dat", virial_history, iterations);
    write_1D_array_to_csv("kinetic_energy_history_T500_P1.dat", kinetic_history, iterations);
    write_1D_array_to_csv("cell_length_history_T500_P1.dat", cell_length_history, iterations);
    //write_2D_array_to_csv("positions_history_20_T500_P1.dat", position_history, iterations, 3);

    destroy_2D_array(positions, n_atoms); 
    destroy_2D_array(velocities, n_atoms);
    destroy_2D_array(forces, n_atoms);    
    destroy_2D_array(position_history, iterations);
}

void simulate_Verlat_task_3(
    double **positions, 
    double **velocities,
    double **forces, 
    double *virials,
    double *potential,
    double mass,
    double isothermal_compression,
    double cell_length,
    double T_eq,
    double P_eq,
    int nbr_atoms,
    int iterations, 
    double time_step,
    double time_constant,
    double *kinetic_history,
    double *virial_history,
    double *cell_length_history, 
    double **positions_history) 
{
    double temp;
    double alpha_t = 1;
    double pressure;
    double alpha_p = 1;

    for (int i = 0; i < iterations; i++) {
        //Calculate the next step.
        velocity_verlet_one_step(
            positions, 
            velocities,
            forces, 
            virials,
            potential,
            mass,
            time_step,
            cell_length,
            nbr_atoms);

        //Save the results from the time step.
        cell_length_history[i] = cell_length;
        virial_history[i] = get_virial_AL(positions, cell_length, nbr_atoms);
        kinetic_history[i] = calculate_kinetic_energy(velocities, mass, nbr_atoms);
        for (int j = 0; j < 3; j++) {
            positions_history[i][j] = positions[20][j];
        }

        //Change the scale of the velocity vector.
        temp = kinetic_history[i] * 2.0 / (3.0 * nbr_atoms * k_B);
        if (i < 3500) {  // Scaling when t < 7 ps.
            alpha_t = sqrt(1 + 2 * time_step/time_constant * (T_eq - temp) / temp);
        } else {         // No scaling when t > 7 ps.
            alpha_t = 1;
        }

        for (int j = 0; j < nbr_atoms; j++) {
            multiplication_with_constant(velocities[j], velocities[j], alpha_t, 3);
        }

        //Change the scale of the position vector.
        pressure = (nbr_atoms * k_B * temp + virial_history[i]) / pow(cell_length, 3);
        if (i < 3500) {   // Scaling when t < 7 ps.
            alpha_p = pow(1 - isothermal_compression * time_step/time_constant * (P_eq - pressure), 1.0/3.0);
        } else {          // No scaling when t > 7 ps.
            alpha_p = 1;
        }

        for (int j = 0; j < nbr_atoms; j++) {
            multiplication_with_constant(positions[j], positions[j], alpha_p, 3);
        }

        //Change the scale of the cell_length
        cell_length *= alpha_p;

        //Print information.
        if ((i+1) % (iterations/20) == 0) {
            printf("Kinetic energy %lf, Temp %lf, alpha_t %lf, ",kinetic_history[i], temp, alpha_t);
            printf("pressure %lf, alpha_p %lf, T_eq %lf, P_eq %lf         |         ", pressure, alpha_p, T_eq, P_eq);
            printf("Progress: %3.1lf %%\n", (double) (i+1)/iterations * 100);
        } 
    }
}

// --------------------------------- Task 4 -------------------------------

void task_4() {
    int n_atoms = 256;
    int N = 4;  // size of supercell in each dimension.
    double mass = 2.7963e-3; // Al mass.
    double lattice_const = 4.;  // Å
    double cell_length = lattice_const * N;

    double bulk_modulus = 70e9 * 6.24e-12; //  eV/Å^3, https://pdf.sciencedirectassets.com/271577/1-s2.0-S0038109800X1455X/1-s2.0-S0038109801005178/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEG4aCXVzLWVhc3QtMSJHMEUCIGeGoeFuBWjXDXfcOs9PKu4Cd9b3xlOaBR%2FrxrrM4Pj%2FAiEA2YNEKeGAywx9m0KMcgc1po0En8SL2Qjz97TIx0nZ3JMqswUIFxAFGgwwNTkwMDM1NDY4NjUiDDVnZR9Zb86Gc5qqtSqQBfrR3NlbxhXrFueJyW6qkO5hTZ7%2BcCTuqPOeSNmFRI26k6aOA4KQzjtZhZi6LnAve3e1a9yN01UYFL%2F4K6nyXL69%2Bsr7%2B5v2pQtxTSq7Q%2BNJXT75D2W64aepYHduL%2BAMDYcLBntpB4D3AxxP3lFDmmmICuOxGILnn2ZWvr0jJPyDkgOcskXRNwriTPiQ%2BbLjJwruGji5HGsDAs7hXLbwDa7EwYVF6ATkyYZMrkOQJ9J027sk8oyAWBX9lHzfXqUm9NsEJDi1kJYEhLrWcHJO42AB613EpiRfO1wSq9rzJF4SmFI%2FE8g3YaYZOal%2B3aH%2BjOhi%2B1q84Bu5fkE%2BWWCUPYnkG8s0aeaTFDFy4O%2FuZ9MHlOHgp3ZOnES8oilMqbdtDM9cNSOqjM4sVB%2Fw5x18qoXNNYWdtex6ix519UEjQUtINd%2FRqOK9I0gpVsId0CIpKvvM5II3RIScmjaJoAMTyCCelxbHWZWWzbR%2BBHroa02MFxekGOiCy6cunDZeNCxAbd%2BXZdi%2BsXnW9lA8G%2F%2FTpn2M6Uac6x6ZtC7rZSzS0Z2%2BcJiQRufh%2BhiSmpVjexUsLXuOpUZbA8XYi32G8%2F2FpK%2FUBepP9jqtAKBUI6YWYkTU8zMxyj3cdTToX8Vn7%2BUc4ut4i6JI9Eq%2Bw8Hh1xkUimu7d5toqbuIDRISspyMGMkF3ml68j1%2FLsNUquT9KV7wi40osXXU0d%2Bs8%2FkLbjeY2jkRsDVrKduTdH40%2FsJlyMBgbV0KowN1SNP2IChjHsUGP%2Bkgb8o7goAHWpzs8MLyVx5jsGLx%2FvxerTTyg6UnhDNJc9Iidvutd4M1WLYDBD1C442HOUVsjYoxzn2JfrnemJ2lBdGvFkfR61bA46xZTsAtMIv4kboGOrEBGADuQ1%2BQlHCvtyEjbIRWPw3xVW30goSR0y0sljN%2BydcVTlnIKaAxN6pcnyzJSyrb%2FHOSDzj9rLx3CY2zXIjIaCBOcc3G7u5gcAfiLS4Vi4ku%2FGqnOaKcazoBglmz4Hj4Nwo4PsdkWyuufrWf5klKVslPYoUIb%2FCjWjI2N6CP7OKmXs0QAC03BQUP5RLJ5HAsyJhqsvR8jldO8FyQQTkChF61zubt%2BSbgT9cddC7tMkN3&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20241125T140307Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYQQGRCNVY%2F20241125%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=ae7fba687d5efc6af8658cbfb805c168673ba6e5061e7d82d0c69fb5aa909772&hash=3ec358775603ee744e95194c31b3e7c635c5ff9545db2af16ba4e34ea442ca9a&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0038109801005178&tid=spdf-314dbf91-09ad-4a81-b820-80c10daef50a&sid=2b46fa2c46952945dc4862b74f1295abb677gxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=15055e015e505d005101&rr=8e822e2bac0febd0&cc=se
    double isothermal_compression = 1/bulk_modulus; // Å^3/eV

    double T_eq_high = 1000. + 273.15;     // Kelvin
    double T_eq_low = 700 + 273.15;
    double P_eq = 1e5 * 6.24e-12;    //eV/Å^3
    
    double **positions = create_2D_array(n_atoms, 3);
    init_fcc(positions, N, lattice_const);
    add_noise_2D_array(positions, 0.065 * lattice_const, n_atoms, 3);

    double **velocities = create_2D_array(n_atoms, 3);
    for (int i = 0; i < n_atoms; i++) {
        for (int j = 0; j < 3; j++) {
            velocities[i][j] = 0;
        }
    }
    double **forces = create_2D_array(n_atoms, 3);

    double virials = 0;
    double potential = 0;

    double time_step = 0.001;
    double time_const = 0.5;
    int iterations_high_temp = 10000;
    int iterations_low_temp = 100000;
    int iterations = iterations_high_temp + iterations_low_temp;

    double kinetic_history[iterations];  
    double virial_history[iterations];
    double cell_length_history[iterations];
    double **position_history = create_2D_array(iterations, 3*n_atoms);
    double **velocities_history = create_2D_array(iterations, 3*n_atoms);
    
     simulate_Verlat_task_4(
        positions, 
        velocities,
        forces, 
        &virials,
        &potential,
        mass,
        isothermal_compression,
        cell_length,
        T_eq_high,
        T_eq_low,
        P_eq,
        n_atoms,
        iterations_high_temp,
        iterations_low_temp, 
        time_step,
        time_const,
        kinetic_history,
        virial_history,
        cell_length_history, 
        position_history,
        velocities_history);

    write_1D_array_to_csv("virial_history_T700_P1_liq.dat", virial_history, iterations);
    write_1D_array_to_csv("kinetic_energy_history_T700_P1_liq.dat", kinetic_history, iterations);
    write_1D_array_to_csv("cell_length_history_T700_P1_liq.dat", cell_length_history, iterations);
    write_2D_array_to_csv("positions_history_T700_P1_liq.dat", position_history, iterations, 3*n_atoms);
    write_2D_array_to_csv("velocities_history_T700_P1_liq.dat", velocities_history, iterations, 3*n_atoms);

    destroy_2D_array(positions, n_atoms); 
    destroy_2D_array(velocities, n_atoms);
    destroy_2D_array(forces, n_atoms);    
    destroy_2D_array(position_history, iterations);
    destroy_2D_array(velocities_history, iterations);
}

void simulate_Verlat_task_4(
    double **positions, 
    double **velocities,
    double **forces, 
    double *virials,
    double *potential,
    double mass,
    double isothermal_compression,
    double cell_length,
    double T_eq_high,
    double T_eq_low,
    double P_eq,
    int nbr_atoms,
    int iterations_high_temp,
    int iterations_low_temp, 
    double time_step,
    double time_constant,
    double *kinetic_history,
    double *virial_history,
    double *cell_length_history, 
    double **positions_history,
    double **velocities_history) 
{
    double temp;
    double alpha_t = 1;
    double pressure;
    double alpha_p = 1;
    int iterations = iterations_high_temp + iterations_low_temp;

    for (int i = 0; i < iterations; i++) {
        //Calculate the next step.
        velocity_verlet_one_step(
            positions, 
            velocities,
            forces, 
            virials,
            potential,
            mass,
            time_step,
            cell_length,
            nbr_atoms);

        //Save the results from the time step.
        cell_length_history[i] = cell_length;
        virial_history[i] = get_virial_AL(positions, cell_length, nbr_atoms);
        kinetic_history[i] = calculate_kinetic_energy(velocities, mass, nbr_atoms);
        for (int k = 0; k < nbr_atoms; ++k) {
            for (int j = 0; j < 3; j++) {
                positions_history[i][j+3*k] = positions[k][j];
            }
        }
        for (int k = 0; k < nbr_atoms; ++k) {
            for (int j = 0; j < 3; j++) {
                velocities_history[i][j+3*k] = velocities[k][j];
            }
        }

        //Change the scale of the velocity vector.
        temp = kinetic_history[i] * 2.0 / (3.0 * nbr_atoms * k_B);
        if (i < iterations_high_temp) {
            alpha_t = sqrt(1 + 2 * time_step/time_constant * (T_eq_high - temp) / temp);
        } else if (i < (iterations_high_temp + 10000)) {
            alpha_t = sqrt(1 + 2 * time_step/time_constant * (T_eq_low - temp) / temp);
        } else {
            alpha_t = 1;
        }
        

        for (int j = 0; j < nbr_atoms; j++) {
            multiplication_with_constant(velocities[j], velocities[j], alpha_t, 3);
        }

        //Change the scale of the position vector.
        pressure = (nbr_atoms * k_B * temp + virial_history[i]) / pow(cell_length, 3);

        if (i < (iterations_high_temp + 10000)) {   // Scaling when t < 7 ps.
            alpha_p = pow(1 - isothermal_compression * time_step/time_constant * (P_eq - pressure), 1.0/3.0);
        } else {                                  // No scaling when t > 7 ps.
            alpha_p = 1;
        }

        for (int j = 0; j < nbr_atoms; j++) {
            multiplication_with_constant(positions[j], positions[j], alpha_p, 3);
        }

        //Change the scale of the cell_length
        cell_length *= alpha_p;

        //Print information.
        if ((i+1) % (iterations/100) == 0) {
            printf("Kinetic energy %lf, Temp %lf, alpha_t %lf, ",kinetic_history[i], temp, alpha_t);
            printf("pressure %lf, alpha_p %lf, T_eq %lf, P_eq %lf         |         ", pressure, alpha_p, T_eq_high, P_eq);
            printf("Progress: %3.1lf %%\n", (double) (i+1)/iterations * 100);
        } 
    }
}
