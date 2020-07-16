// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

inline void indirect_rw_kernel(
    #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
        double ewt,
    #endif
    double ex, double ey, double ez,
    #if defined SIMD && defined MANUAL_GATHER
        const double variables_a[][DBLS_PER_SIMD],
        const double variables_b[][DBLS_PER_SIMD],
    #else
        const double *restrict variables_a, 
        const double *restrict variables_b, 
    #endif
    #ifdef FLUX_FISSION
        edge *restrict edge_variables
    #else
        #if defined SIMD && defined MANUAL_SCATTER
            int simd_idx,
            double fluxes_a[][DBLS_PER_SIMD],
            double fluxes_b[][DBLS_PER_SIMD]
        #else
            double *restrict fluxes_a, 
            double *restrict fluxes_b
        #endif
    #endif
    )
{
    // Process edge-point A:
    double p_a, pe_a;
    double3 momentum_a;
    const int p_a_idx  = VAR_DENSITY;
    const int mx_a_idx = VAR_MOMENTUMX;
    const int my_a_idx = VAR_MOMENTUMY;
    const int mz_a_idx = VAR_MOMENTUMZ;
    const int pe_a_idx = VAR_DENSITY_ENERGY;
    #if defined SIMD && defined MANUAL_GATHER
        p_a          = variables_a[ p_a_idx][simd_idx];
        momentum_a.x = variables_a[mx_a_idx][simd_idx];
        momentum_a.y = variables_a[my_a_idx][simd_idx];
        momentum_a.z = variables_a[mz_a_idx][simd_idx];
        pe_a         = variables_a[pe_a_idx][simd_idx];
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
    const int p_b_idx  = VAR_DENSITY;
    const int mx_b_idx = VAR_MOMENTUMX;
    const int my_b_idx = VAR_MOMENTUMY;
    const int mz_b_idx = VAR_MOMENTUMZ;
    const int pe_b_idx = VAR_DENSITY_ENERGY;
    #if defined SIMD && defined MANUAL_GATHER
        p_b          = variables_b[ p_b_idx][simd_idx];
        momentum_b.x = variables_b[mx_b_idx][simd_idx];
        momentum_b.y = variables_b[my_b_idx][simd_idx];
        momentum_b.z = variables_b[mz_b_idx][simd_idx];
        pe_b         = variables_b[pe_b_idx][simd_idx];
    #else
        p_b          = variables_b[ p_b_idx];
        momentum_b.x = variables_b[mx_b_idx];
        momentum_b.y = variables_b[my_b_idx];
        momentum_b.z = variables_b[mz_b_idx];
        pe_b         = variables_b[pe_b_idx];
    #endif

    double p_a_val  = p_b + ex;
    double pe_a_val = pe_b + ey;
    double mx_a_val = momentum_b.x + ez;
    double my_a_val = momentum_b.y;
    double mz_a_val = momentum_b.z;
    #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
        double p_b_val = p_a + ewt;
        #else
        double p_b_val = p_a;
    #endif
    double pe_b_val = pe_a;
    double mx_b_val = momentum_a.x;
    double my_b_val = momentum_a.y;
    double mz_b_val = momentum_a.z;

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
            fluxes_a[VAR_DENSITY][simd_idx]   += p_a_val;
            fluxes_a[VAR_MOMENTUMX][simd_idx] += mx_a_val;
            fluxes_a[VAR_MOMENTUMY][simd_idx] += my_a_val;
            fluxes_a[VAR_MOMENTUMZ][simd_idx] += mz_a_val;
            fluxes_a[VAR_DENSITY_ENERGY][simd_idx] += pe_a_val;

            fluxes_b[VAR_DENSITY][simd_idx]   += p_b_val;
            fluxes_b[VAR_MOMENTUMX][simd_idx] += mx_b_val;
            fluxes_b[VAR_MOMENTUMY][simd_idx] += my_b_val;
            fluxes_b[VAR_MOMENTUMZ][simd_idx] += mz_b_val;
            fluxes_b[VAR_DENSITY_ENERGY][simd_idx] += pe_b_val;
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
