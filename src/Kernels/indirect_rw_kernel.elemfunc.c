// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

/*
inline void indirect_rw_kernel(
    #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
        double ewt,
    #endif
    #if defined SIMD && defined MANUAL_GATHER
        const double edge_vector[][DBLS_PER_SIMD],
        const double variables_a[][DBLS_PER_SIMD],
        const double variables_b[][DBLS_PER_SIMD],
    #else
        double ex, double ey, double ez,
        const double *restrict variables_a, 
        const double *restrict variables_b, 
    #endif
    #ifdef FLUX_FISSION
        edge *restrict edge_variables
    #else
        #if defined SIMD && defined MANUAL_SCATTER
            int simd_idx,
            double fluxes_a[][MANUAL_WIDTH],
            double fluxes_b[][MANUAL_WIDTH]
        #else
            double *restrict fluxes_a, 
            double *restrict fluxes_b
        #endif
    #endif
    )
*/

    const long a = edges[i].a;
    const long b = edges[i].b;
    #if defined SIMD && (defined MANUAL_GATHER || defined MANUAL_SCATTER)
        const int simd_idx = i - loop_start;
    #endif

    #if defined SIMD && (defined MANUAL_GATHER || defined MANUAL_SCATTER)
        double ex = edge_vectors[0][simd_idx];
        double ey = edge_vectors[1][simd_idx];
        double ez = edge_vectors[2][simd_idx];
    #else
        double ex = edges[i].x;
        double ey = edges[i].y;
        double ez = edges[i].z;
    #endif
    #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
        double ewt = edge_weights[i];
    #endif

    ////////////////////////////////////
    // Read and process edge-point B:
    ////////////////////////////////////
    double p_b, pe_b;
    double3 momentum_b;
    const long p_b_idx  = b*NVAR + VAR_DENSITY;
    const long mx_b_idx = b*NVAR + VAR_MOMENTUMX;
    const long my_b_idx = b*NVAR + VAR_MOMENTUMY;
    const long mz_b_idx = b*NVAR + VAR_MOMENTUMZ;
    const long pe_b_idx = b*NVAR + VAR_DENSITY_ENERGY;
    #if defined SIMD && defined MANUAL_GATHER
        p_b          = variables_b[VAR_DENSITY]       [simd_idx];
        momentum_b.x = variables_b[VAR_MOMENTUMX]     [simd_idx];
        momentum_b.y = variables_b[VAR_MOMENTUMY]     [simd_idx];
        momentum_b.z = variables_b[VAR_MOMENTUMZ]     [simd_idx];
        pe_b         = variables_b[VAR_DENSITY_ENERGY][simd_idx];
    #else
        p_b          = variables[ p_b_idx];
        momentum_b.x = variables[mx_b_idx];
        momentum_b.y = variables[my_b_idx];
        momentum_b.z = variables[mz_b_idx];
        pe_b         = variables[pe_b_idx];
    #endif

    ////////////////////////////////////
    // Read and process edge-point A:
    ////////////////////////////////////
    double p_a, pe_a;
    double3 momentum_a;
    const long p_a_idx  = a*NVAR + VAR_DENSITY;
    const long mx_a_idx = a*NVAR + VAR_MOMENTUMX;
    const long my_a_idx = a*NVAR + VAR_MOMENTUMY;
    const long mz_a_idx = a*NVAR + VAR_MOMENTUMZ;
    const long pe_a_idx = a*NVAR + VAR_DENSITY_ENERGY;
    #if defined SIMD && defined MANUAL_GATHER
        p_a          = variables_a[VAR_DENSITY]       [simd_idx];
        momentum_a.x = variables_a[VAR_MOMENTUMX]     [simd_idx];
        momentum_a.y = variables_a[VAR_MOMENTUMY]     [simd_idx];
        momentum_a.z = variables_a[VAR_MOMENTUMZ]     [simd_idx];
        pe_a         = variables_a[VAR_DENSITY_ENERGY][simd_idx];
    #else
        p_a          = variables[ p_a_idx];
        momentum_a.x = variables[mx_a_idx];
        momentum_a.y = variables[my_a_idx];
        momentum_a.z = variables[mz_a_idx];
        pe_a         = variables[pe_a_idx];
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
        edge_variables[i*NVAR + VAR_DENSITY       ].a =  p_a_val;
        edge_variables[i*NVAR + VAR_MOMENTUMX     ].a = mx_a_val;
        edge_variables[i*NVAR + VAR_MOMENTUMY     ].a = my_a_val;
        edge_variables[i*NVAR + VAR_MOMENTUMZ     ].a = mz_a_val;
        edge_variables[i*NVAR + VAR_DENSITY_ENERGY].a = pe_a_val;

        edge_variables[i*NVAR + VAR_DENSITY       ].b =  p_b_val;
        edge_variables[i*NVAR + VAR_MOMENTUMX     ].b = mx_b_val;
        edge_variables[i*NVAR + VAR_MOMENTUMY     ].b = my_b_val;
        edge_variables[i*NVAR + VAR_MOMENTUMZ     ].b = mz_b_val;
        edge_variables[i*NVAR + VAR_DENSITY_ENERGY].b = pe_b_val;
    #else
        const long p_a_flx_idx  = a*NVAR + VAR_DENSITY;
        const long mx_a_flx_idx = a*NVAR + VAR_MOMENTUMX;
        const long my_a_flx_idx = a*NVAR + VAR_MOMENTUMY;
        const long mz_a_flx_idx = a*NVAR + VAR_MOMENTUMZ;
        const long pe_a_flx_idx = a*NVAR + VAR_DENSITY_ENERGY;

        const long p_b_flx_idx  = b*NVAR + VAR_DENSITY;
        const long mx_b_flx_idx = b*NVAR + VAR_MOMENTUMX;
        const long my_b_flx_idx = b*NVAR + VAR_MOMENTUMY;
        const long mz_b_flx_idx = b*NVAR + VAR_MOMENTUMZ;
        const long pe_b_flx_idx = b*NVAR + VAR_DENSITY_ENERGY;

        #if defined USE_AVX512CD
            #pragma omp ordered simd overlap(p_a_flx_idx)
            {
                fluxes[p_a_flx_idx]  +=  p_a_val;
            }
            #pragma omp ordered simd overlap(mx_a_flx_idx)
            {
                fluxes[mx_a_flx_idx] += mx_a_val;
            }
            #pragma omp ordered simd overlap(my_a_flx_idx)
            {
                fluxes[my_a_flx_idx] += my_a_val;
            }
            #pragma omp ordered simd overlap(mz_a_flx_idx)
            {
                fluxes[mz_a_flx_idx] += mz_a_val;
            }
            #pragma omp ordered simd overlap(pe_a_flx_idx)
            {
                fluxes[pe_a_flx_idx] += pe_a_val;
            }

            #pragma omp ordered simd overlap(p_b_flx_idx)
            {
                fluxes[p_b_flx_idx]  +=  p_b_val;
            }
            #pragma omp ordered simd overlap(mx_b_flx_idx)
            {
                fluxes[mx_b_flx_idx] += mx_b_val;
            }
            #pragma omp ordered simd overlap(my_b_flx_idx)
            {
                fluxes[my_b_flx_idx] += my_b_val;
            }
            #pragma omp ordered simd overlap(mz_b_flx_idx)
            {
                fluxes[mz_b_flx_idx] += mz_b_val;
            }
            #pragma omp ordered simd overlap(pe_b_flx_idx)
            {
                fluxes[pe_b_flx_idx] += pe_b_val;
            }
        #else
            #if defined SIMD && defined MANUAL_SCATTER
                fluxes_a[VAR_DENSITY]  [simd_idx]      = p_a_val;
                fluxes_a[VAR_MOMENTUMX][simd_idx]      = mx_a_val;
                fluxes_a[VAR_MOMENTUMY][simd_idx]      = my_a_val;
                fluxes_a[VAR_MOMENTUMZ][simd_idx]      = mz_a_val;
                fluxes_a[VAR_DENSITY_ENERGY][simd_idx] = pe_a_val;

                fluxes_b[VAR_DENSITY]  [simd_idx]      = p_b_val;
                fluxes_b[VAR_MOMENTUMX][simd_idx]      = mx_b_val;
                fluxes_b[VAR_MOMENTUMY][simd_idx]      = my_b_val;
                fluxes_b[VAR_MOMENTUMZ][simd_idx]      = mz_b_val;
                fluxes_b[VAR_DENSITY_ENERGY][simd_idx] = pe_b_val;
            #else
                fluxes[p_a_flx_idx]  +=  p_a_val;
                fluxes[mx_a_flx_idx] += mx_a_val;
                fluxes[my_a_flx_idx] += my_a_val;
                fluxes[mz_a_flx_idx] += mz_a_val;
                fluxes[pe_a_flx_idx] += pe_a_val;

                fluxes[p_b_flx_idx]  +=  p_b_val;
                fluxes[mx_b_flx_idx] += mx_b_val;
                fluxes[my_b_flx_idx] += my_b_val;
                fluxes[mz_b_flx_idx] += mz_b_val;
                fluxes[pe_b_flx_idx] += pe_b_val;
            #endif
        #endif
    #endif
