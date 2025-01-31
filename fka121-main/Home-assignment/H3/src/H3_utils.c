#include "H3_utils.h"

double
pot_task1(double x)
{
    return 1./2. * pow(1. - exp(-x), 2.);
}

double 
pot_task2(double *x)
{
    double pos1[] = {x[0], x[1], x[2]};
    double pos2[] = {x[3], x[4], x[5]};

    return -2 * (1/vector_norm(pos1, 3) + 1/vector_norm(pos2, 3)) + 1/distance_between_vectors(pos1, pos2, 3); 
}
 
double calculate_average_E_T(
    int numbSamples,
    double new_E_T,
    double average_E_T) 
{
    return 1. / (numbSamples + 1.) * new_E_T + numbSamples / (1. + numbSamples) * average_E_T;
}

double calculate_W(
    double x,
    double delta_t,
    double E_T, 
    double (*potential)(double)
)
{
    double W = exp(-((*potential)(x) - E_T) * delta_t);
    return W;
}

double calculate_W2(
    double *x,
    double delta_t,
    double E_T, 
    double (*potential)(double *)
)
{
    double W = exp(-((*potential)(x) - E_T) * delta_t);
    return W;
}

int numberOfWalkers(
    double x,
    double delta_t,
    double E_T, 
    double (*potential)(double),
    gsl_rng *r) 
{
    double W = exp(-((*potential)(x) - E_T) * delta_t);
    double U = gsl_rng_uniform(r);
    return (int)(W + U);
}

int numberOfWalkers2(
    double *x,
    double delta_t,
    double E_T, 
    double (*potential)(double *),
    gsl_rng *r) 
{
    double W = exp(-((*potential)(x) - E_T) * delta_t);
    double U = gsl_rng_uniform(r);
    return (int)(W + U);
}
