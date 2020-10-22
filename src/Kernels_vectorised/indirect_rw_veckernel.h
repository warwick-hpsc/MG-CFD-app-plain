// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

#ifndef INDIRECT_RW_VECKERNEL_H
#define INDIRECT_RW_VECKERNEL_H

FORCE_INLINE
inline void indirect_rw_veckernel(
    #if defined SIMD && (defined MANUAL_GATHER || defined MANUAL_SCATTER)
        int simd_idx,
    #endif

    #if defined SIMD && defined MANUAL_GATHER
        #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
            const double simd_edge_weights[DBLS_PER_SIMD],
        #else
        #endif
        const double simd_edge_vectors[][DBLS_PER_SIMD],
        const double simd_variables_a[][DBLS_PER_SIMD],
        const double simd_variables_b[][DBLS_PER_SIMD],
    #else
        #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
            double ewt,
        #endif
        double ex, double ey, double ez,
        const double *restrict variables_a, 
        const double *restrict variables_b, 
    #endif

    #if defined SIMD && defined MANUAL_SCATTER
        double simd_fluxes_a[][DBLS_PER_SIMD],
        double simd_fluxes_b[][DBLS_PER_SIMD]
    #elif defined FLUX_FISSION
        edge *restrict edge_variables
    #else
        double *restrict fluxes_a, 
        double *restrict fluxes_b
    #endif
    )
{
    #if defined SIMD && defined MANUAL_GATHER
        double ex = simd_edge_vectors[0][simd_idx];
        double ey = simd_edge_vectors[1][simd_idx];
        double ez = simd_edge_vectors[2][simd_idx];
        #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
            double ewt = edge_weights[simd_idx];
        #else
            double ewt = sqrt(ex*ex + ey*ey + ez*ez);
        #endif
    #else
        #ifndef FLUX_PRECOMPUTE_EDGE_WEIGHTS
            double ewt = sqrt(ex*ex + ey*ey + ez*ez);
        #endif
    #endif

    // Process edge-point A:
    double p_a, pe_a;
    double3 momentum_a;
    const long p_a_idx  = VAR_DENSITY;
    const long mx_a_idx = VAR_MOMENTUMX;
    const long my_a_idx = VAR_MOMENTUMY;
    const long mz_a_idx = VAR_MOMENTUMZ;
    const long pe_a_idx = VAR_DENSITY_ENERGY;
    #if defined SIMD && defined MANUAL_GATHER
        p_a          = simd_variables_a[ p_a_idx][simd_idx];
        momentum_a.x = simd_variables_a[mx_a_idx][simd_idx];
        momentum_a.y = simd_variables_a[my_a_idx][simd_idx];
        momentum_a.z = simd_variables_a[mz_a_idx][simd_idx];
        pe_a         = simd_variables_a[pe_a_idx][simd_idx];
    #else
        p_a          = variables_a[ p_a_idx];
        momentum_a.x = variables_a[mx_a_idx];
        momentum_a.y = variables_a[my_a_idx];
        momentum_a.z = variables_a[mz_a_idx];
        pe_a         = variables_a[pe_a_idx];
    #endif

    // Process edge-point B:
    double p_b, pe_b;
    double3 momentum_b;
    const long p_b_idx  = VAR_DENSITY;
    const long mx_b_idx = VAR_MOMENTUMX;
    const long my_b_idx = VAR_MOMENTUMY;
    const long mz_b_idx = VAR_MOMENTUMZ;
    const long pe_b_idx = VAR_DENSITY_ENERGY;
    #if defined SIMD && defined MANUAL_GATHER
        p_b          = simd_variables_b[ p_b_idx][simd_idx];
        momentum_b.x = simd_variables_b[mx_b_idx][simd_idx];
        momentum_b.y = simd_variables_b[my_b_idx][simd_idx];
        momentum_b.z = simd_variables_b[mz_b_idx][simd_idx];
        pe_b         = simd_variables_b[pe_b_idx][simd_idx];
    #else
        p_b          = variables_b[ p_b_idx];
        momentum_b.x = variables_b[mx_b_idx];
        momentum_b.y = variables_b[my_b_idx];
        momentum_b.z = variables_b[mz_b_idx];
        pe_b         = variables_b[pe_b_idx];
    #endif

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
        #if defined SIMD && defined MANUAL_SCATTER
            simd_fluxes_a[VAR_DENSITY][simd_idx]        = p_a_val;
            simd_fluxes_a[VAR_MOMENTUMX][simd_idx]      = mx_a_val;
            simd_fluxes_a[VAR_MOMENTUMY][simd_idx]      = my_a_val;
            simd_fluxes_a[VAR_MOMENTUMZ][simd_idx]      = mz_a_val;
            simd_fluxes_a[VAR_DENSITY_ENERGY][simd_idx] = pe_a_val;

            simd_fluxes_b[VAR_DENSITY][simd_idx]        = p_b_val;
            simd_fluxes_b[VAR_MOMENTUMX][simd_idx]      = mx_b_val;
            simd_fluxes_b[VAR_MOMENTUMY][simd_idx]      = my_b_val;
            simd_fluxes_b[VAR_MOMENTUMZ][simd_idx]      = mz_b_val;
            simd_fluxes_b[VAR_DENSITY_ENERGY][simd_idx] = pe_b_val;
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
    #endif
}

#endif
