#ifndef VALIDATION_H
#define VALIDATION_H

#include "definitions.h"

void adjust_ewt(
    const double3 *restrict coords,
    long num_edges, 
    edge_neighbour *restrict edges);

void dampen_ewt(
    long num_edges, 
    edge_neighbour *restrict edges, 
    double damping_factor);

void residual(
    long nel, 
    const double *restrict old_variables, 
    const double *restrict variables, 
    double *restrict residuals);

double calc_rms(
    long nel, 
    const double *restrict residuals);

void check_for_invalid_variables(
    const double *restrict variables,
    long n);

void identify_differences(
    const double* test_values,
    const double* master_values, 
    long n, 
    const long* peritab);

#endif
