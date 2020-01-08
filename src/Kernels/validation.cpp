//************************************************//
// Copyright 2016-2019 University of Warwick

// Permission is hereby granted, free of charge, to any person obtaining 
// a copy of this software and associated documentation files (the "Software"), 
// to deal in the Software without restriction, including without limitation 
// the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or 
// sell copies of the Software, and to permit persons to whom the Software is furnished 
// to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included 
// in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
// PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
// COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
// IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//************************************************//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "validation.h"
#include "common.h"

void adjust_ewt(
    const double3 *restrict coords, 
    long num_edges, 
    edge_neighbour *restrict edges)
{
    // |ewt| currently is face area. Divide through by distance 
    // to produce 'surface vector' with magnitude (area/dm):

    log("adjust_ewt()\n");
    #ifdef OMP
        #pragma omp parallel for
    #endif
    for (long i=0; i<num_edges; i++) {
        long a = edges[i].a;
        long b = edges[i].b;

        if (a >= 0 && b >= 0) {
            double dist = 0.0, d;
            d = coords[b].x - coords[a].x;
            dist += d*d;
            d = coords[b].y - coords[a].y;
            dist += d*d;
            d = coords[b].z - coords[a].z;
            dist += d*d;
            dist = sqrt(dist);

            edges[i].x /= dist;
            edges[i].y /= dist;
            edges[i].z /= dist;
        }
    }
}

void dampen_ewt(
    long num_edges, 
    edge_neighbour *restrict edges, 
    double damping_factor)
{
    log("dampen_ewt()\n");
    #ifdef OMP
        #pragma omp parallel for
    #endif
    for (long i=0; i<num_edges; i++) {
        edges[i].x *= damping_factor;
        edges[i].y *= damping_factor;
        edges[i].z *= damping_factor;
    }
}

void residual(
    long nel, 
    const double *restrict old_variables, 
    const double *restrict variables, 
    double *restrict residuals)
{
    #ifdef OMP
        #pragma omp parallel for
    #endif
    for (long i=0; i<(nel*NVAR); i++) {
        residuals[i] = variables[i] - old_variables[i];
    }
}

double calc_rms(
    long nel, 
    const double *restrict residuals)
{
    double rms = 0.0;
    #ifdef OMP
        #pragma omp parallel for reduction(+:rms)
    #endif
    for (long i=0; i<(nel*NVAR); i++) {
        rms += pow(residuals[i], 2);
    }
    rms /= double(nel);
    rms = sqrt(rms);
    return rms;
}

void check_for_invalid_variables(
    const double *restrict variables,
    long n)
{   
    for (long i=0; i<n; i++) {
        for (int v=0; v<NVAR; v++) {
            const long idx = i*NVAR + v;

            if (isnan(variables[idx]) || isinf(variables[idx])) {
                printf("\nERROR: NaN detected!");
                printf("\nCell %ld: d  = %.4e", i, variables[i*NVAR + 0]);
                printf("\nCell %ld: mx = %.4e", i, variables[i*NVAR + 1]);
                printf("\nCell %ld: my = %.4e", i, variables[i*NVAR + 2]);
                printf("\nCell %ld: mz = %.4e", i, variables[i*NVAR + 3]);
                printf("\nCell %ld: de = %.4e", i, variables[i*NVAR + 4]);
                exit(EXIT_FAILURE);
            }
        }

        if (variables[i*NVAR + VAR_DENSITY] < 0.0) {
            printf("\nERROR: Negative density detected!");
            printf("\nCell %ld: d  = %.4e", i, variables[i*NVAR + VAR_DENSITY]);
            exit(EXIT_FAILURE);
        }

        if (variables[i*NVAR + VAR_DENSITY_ENERGY] < 0.0) {
            printf("\nERROR: Negative density.energy detected!");
            printf("\nCell %ld: de  = %.4e", i, variables[i*NVAR + VAR_DENSITY_ENERGY]);
            exit(EXIT_FAILURE);
        }
    }
}

void identify_differences(
    const double* test_values,
    const double* master_values, 
    long n)
{
    // If floating-point operations have been reordered, then a difference
    // is expected due to rounding-errors, but the difference should
    // be smaller than the following:
    //   1 x 10 ^ ( E - 17 + N )
    // Where E = exponent of master value
    //       N = largest difference in exponents of any floating-point
    //           arithmetic operation performed

    // N represents how many of the least-significant base-10 digits
    // of the floating-point mantissa are allowed to differ due to 
    // FP arithmetic reordering. Its value is guessed as 8, as to set it
    // accurately would require a trace of all floating-point operation
    // outputs during the runs.

    const double acceptable_relative_difference = 10.0e-9;

    double absolute_threshold = 3.0e-19;
    if (mesh_variant == MESH_FVCORR) {
        // Relax threshold for this mesh, as original code performs 
        // arithmetic in a hugely different order.
        absolute_threshold = 1.0e-15;
    }

    for (long i=0; i<n; i++) {
        for (int v=0; v<NVAR; v++) {
            const long idx = i*NVAR + v;

            double acceptable_difference = master_values[idx] * acceptable_relative_difference;
            if (acceptable_difference < 0.0) {
                acceptable_difference *= -1.0;
            }

            // Ignore any differences smaller than 3e-19:
            if (acceptable_difference < absolute_threshold) {
                acceptable_difference = absolute_threshold;
            }

            double diff = test_values[idx] - master_values[idx];
            if (diff < 0.0) {
                diff *= -1.0;
            }

            if (diff > acceptable_difference) {
                printf("ERROR: Unacceptable error detected at (i=%ld, v=%d)\n", i, v);
                // printf("       - incorrect value = %.17f\n", test_values[idx]);
                // printf("       - correct value =   %.17f\n", master_values[idx]);
                // printf("       - diff          =   %.17f\n", diff);
                printf("       - incorrect value = %.23f\n", test_values[idx]);
                printf("       - correct value =   %.23f\n", master_values[idx]);
                printf("       - diff          =   %.23f\n", diff);
                exit(EXIT_FAILURE);
            }
        }
    }
}
