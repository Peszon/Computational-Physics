#pragma once
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tools.h"

/* **********************************************
 *
 * Returns the value of the 1-D Morse potential.
 *
 * **********************************************/
double
pot_task1(double x);


/* **********************************************
 *
 * Returns the value of the 3-D Helium potential.
 *
 * **********************************************/
double
pot_task2(double *x);

/* **********************************************
 *
 * Returns the new value of the average E_T.
 *
 * **********************************************/
double calculate_average_E_T(
    int numbSamples,
    double new_E_T,
    double average_E_T); 


/* **********************************************
 *
 * Returns the new value of W of the walker with 
 * the current average energy and the walkers 
 * position, task 1.
 *
 * **********************************************/
double calculate_W(
    double x,
    double delta_t,
    double E_T, 
    double (*potential)(double));


/* **********************************************
 *
 * Returns the new value of W of the walker with 
 * the current average energy and the walkers 
 * position, task 2.
 *
 * **********************************************/
double calculate_W2(
    double *x,
    double delta_t,
    double E_T, 
    double (*potential)(double *)
);

/* **********************************************
 * Task 1.
 * Returns the amount of walkers in the next timestep.
 * 
 * 0 = decrese, 1 = no change, 2,3... = Increase.
 *
 * **********************************************/
int numberOfWalkers(
    double x,
    double delta_t,
    double E_T, 
    double (*potential)(double),
    gsl_rng *r);


/* **********************************************
 * Task 2.
 * Returns the amount of walkers in the next timestep.
 * 
 * 0 = decrese, 1 = no change, 2,3... = Increase.
 *
 * **********************************************/
int numberOfWalkers2(
    double *x,
    double delta_t,
    double E_T, 
    double (*potential)(double *),
    gsl_rng *r);
