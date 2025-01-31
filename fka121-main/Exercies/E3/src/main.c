#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "tools.h"

#define PI 3.141592653589
#define k_B 1.380649e-8  // μm² μg/ms² K

typedef struct {
    double position;
    double velocity; 
} result_t;

void task3();
result_t BD3(double initial_position, double initial_velocity, double w0, double dt, double eta, double kB, double mass, double T, gsl_rng *k);

int main() {
    task3();
    return 1;
}

void task3(){
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    double w0 = 3.1 * 2 * PI;
    double dt[] = {0.001, 0.005}; 
    double eta[] = {1/0.1473, 1/0.0485};  // milli second^-1 
    double mass = 3.01305e-5;   // micro g 

    int n_timesteps = 0.1473 / dt[0] * 5000;

    double **position_histories = create_2D_array(n_timesteps, 4);
    double **velocity_histories = create_2D_array(n_timesteps, 4);

    for(int i = 0; i < 2; i++) {
        double relaxation_time = 1/eta[i]; 
        for(int j = 0; j < 2; j++) {
            position_histories[0][i*2 + j] = 0.04;
            velocity_histories[0][i*2 + j] = 0.1;
            for(int k = 1; k < n_timesteps; k++) {
                result_t results = BD3(position_histories[k-1][i*2 + j], velocity_histories[k-1][i*2 + j], w0, dt[j], eta[i], k_B, mass, 300, r);
                position_histories[k][i*2 + j] = results.position;
                velocity_histories[k][i*2 + j] = results.velocity;
            }
        }
    }

    printf("Number of timesteps: %d\n", n_timesteps);

    write_2D_array_to_csv("results/task3/position_data.dat", position_histories, n_timesteps, 4);
    write_2D_array_to_csv("results/task3/velocities_data.dat", velocity_histories, n_timesteps, 4);

    destroy_2D_array(position_histories, n_timesteps);
    destroy_2D_array(velocity_histories, n_timesteps);
    gsl_rng_free(r);
}

result_t BD3(double initial_position, double initial_velocity, double w0, double dt, double eta, double kB, double mass, double T, gsl_rng *k)
{
    result_t result;

    double c0 = exp(-eta * dt);
    double sqrt_c0 = pow(c0, 1./2.);
    double vth = 0.3707653151;
 
    double G1 = gsl_ran_gaussian(k, 1.0);
    double G2 = gsl_ran_gaussian(k, 1.0);

    double a_n = -pow(w0, 2) * initial_position;

    double vtilde_n = 1./2. * a_n * dt + sqrt_c0 * initial_velocity + pow(1-c0, 1./2.) * vth * G1;

    result.position = initial_position + vtilde_n * dt;

    double a_n2 = -pow(w0, 2) * result.position;

    result.velocity = sqrt_c0 * a_n2 * dt/2 + sqrt_c0 * vtilde_n + vth * pow(1-c0, 1./2.) * G2;

    return result;
}