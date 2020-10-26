// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

#ifndef UNSTRUCTURED_STREAM_KERNEL_H
#define UNSTRUCTURED_STREAM_KERNEL_H

FORCE_INLINE
inline void unstructured_stream_kernel(
    #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
        double ewt,
    #endif
    double ex, double ey, double ez,
    const double *restrict variables_a, 
    const double *restrict variables_b, 
    #ifdef FLUX_FISSION
        edge *restrict edge_variables
    #else
        double *restrict fluxes_a, 
        double *restrict fluxes_b
    #endif
    )
{
    #ifndef FLUX_PRECOMPUTE_EDGE_WEIGHTS
        double ewt = sqrt(ex*ex + ey*ey + ez*ez);
    #endif

    // Process edge-point A:
    double p_a, pe_a;
    double3 momentum_a;
    const long p_a_idx  = VAR_DENSITY;
    const long mx_a_idx = VAR_MOMENTUMX;
    const long my_a_idx = VAR_MOMENTUMY;
    const long mz_a_idx = VAR_MOMENTUMZ;
    const long pe_a_idx = VAR_DENSITY_ENERGY;
    p_a          = variables_a[ p_a_idx];
    momentum_a.x = variables_a[mx_a_idx];
    momentum_a.y = variables_a[my_a_idx];
    momentum_a.z = variables_a[mz_a_idx];
    pe_a         = variables_a[pe_a_idx];

    // Process edge-point B:
    double p_b, pe_b;
    double3 momentum_b;
    const long p_b_idx  = VAR_DENSITY;
    const long mx_b_idx = VAR_MOMENTUMX;
    const long my_b_idx = VAR_MOMENTUMY;
    const long mz_b_idx = VAR_MOMENTUMZ;
    const long pe_b_idx = VAR_DENSITY_ENERGY;
    p_b          = variables_b[ p_b_idx];
    momentum_b.x = variables_b[mx_b_idx];
    momentum_b.y = variables_b[my_b_idx];
    momentum_b.z = variables_b[mz_b_idx];
    pe_b         = variables_b[pe_b_idx];

    double p_a_val  = p_b + ex;
    double mx_a_val = momentum_b.y * momentum_b.z + ez;
    double my_a_val = momentum_b.x * momentum_b.z;
    double mz_a_val = momentum_b.x + momentum_b.y;
    double pe_a_val = pe_b - ey;
    
    #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
        double p_b_val = p_a - ewt;
        #else
        double p_b_val = p_a - ey;
    #endif
    double mx_b_val = momentum_a.y + momentum_a.z;
    double my_b_val = momentum_a.x * momentum_a.z - ez;
    double mz_b_val = momentum_a.x - momentum_a.y;
    double pe_b_val = pe_a + ey;
    
    // Write out fluxes to memory:
    #ifdef FLUX_FISSION
        edge_variables[VAR_DENSITY       ].a =  p_a_val;
        edge_variables[VAR_MOMENTUMX     ].a = mx_a_val;
        edge_variables[VAR_MOMENTUMY     ].a = my_a_val;
        edge_variables[VAR_MOMENTUMZ     ].a = mz_a_val;
        edge_variables[VAR_DENSITY_ENERGY].a = pe_a_val;

        edge_variables[VAR_DENSITY       ].b =  p_b_val;
        edge_variables[VAR_MOMENTUMX     ].b = mx_b_val;
        edge_variables[VAR_MOMENTUMY     ].b = my_b_val;
        edge_variables[VAR_MOMENTUMZ     ].b = mz_b_val;
        edge_variables[VAR_DENSITY_ENERGY].b = pe_b_val;
    #else
        fluxes_a[VAR_DENSITY]  +=  p_a_val;
        fluxes_a[VAR_MOMENTUMX] += mx_a_val;
        fluxes_a[VAR_MOMENTUMY] += my_a_val;
        fluxes_a[VAR_MOMENTUMZ] += mz_a_val;
        fluxes_a[VAR_DENSITY_ENERGY] += pe_a_val;

        fluxes_b[VAR_DENSITY]  +=  p_b_val;
        fluxes_b[VAR_MOMENTUMX] += mx_b_val;
        fluxes_b[VAR_MOMENTUMY] += my_b_val;
        fluxes_b[VAR_MOMENTUMZ] += mz_b_val;
        fluxes_b[VAR_DENSITY_ENERGY] += pe_b_val;
    #endif
}

#endif
