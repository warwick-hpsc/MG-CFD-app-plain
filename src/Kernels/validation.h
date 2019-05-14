#ifndef VALIDATION_H
#define VALIDATION_H

#include "definitions.h"

void residual(
    int nel, 
    const double *restrict old_variables, 
    const double *restrict variables, 
    double *restrict residuals);

double calc_rms(
    int nel, 
    const double *restrict residuals);

void check_for_invalid_variables(
    const double *restrict variables,
    int n);

void identify_differences(
    const double* test_values,
    const double* master_values, 
    int n);

#endif
