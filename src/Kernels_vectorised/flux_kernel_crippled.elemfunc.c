// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

/*
inline void compute_flux_edge_kernel_crippled(
    #if defined SIMD && (defined MANUAL_GATHER || defined MANUAL_SCATTER)
        int simd_idx,
    #endif

    #if defined SIMD && defined MANUAL_GATHER
        #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
            const double simd_edge_weights[DBLS_PER_SIMD],
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
*/

    const long a = edge_nodes[i*2];
    const long b = edge_nodes[i*2+1];
    #if defined SIMD && (defined MANUAL_GATHER || defined MANUAL_SCATTER)
        const int simd_idx = i - flux_loop_start;
    #endif

    #if defined SIMD && defined MANUAL_GATHER
        double ex = simd_edge_vectors[0][simd_idx];
        double ey = simd_edge_vectors[1][simd_idx];
        double ez = simd_edge_vectors[2][simd_idx];
        #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
            double ewt = simd_edge_weights[0][simd_idx];
        #else
            double ewt = sqrt(ex*ex + ey*ey + ez*ez);
        #endif
    #else
        double ex = edge_vectors[i*NDIM];
        double ey = edge_vectors[i*NDIM+1];
        double ez = edge_vectors[i*NDIM+2];
        #ifndef FLUX_PRECOMPUTE_EDGE_WEIGHTS
            double ewt = sqrt(ex*ex + ey*ey + ez*ez);
        #endif
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
        p_b          = simd_variables_b[VAR_DENSITY]       [simd_idx];
        momentum_b.x = simd_variables_b[VAR_MOMENTUMX]     [simd_idx];
        momentum_b.y = simd_variables_b[VAR_MOMENTUMY]     [simd_idx];
        momentum_b.z = simd_variables_b[VAR_MOMENTUMZ]     [simd_idx];
        pe_b         = simd_variables_b[VAR_DENSITY_ENERGY][simd_idx];
    #else
        p_b          = variables[ p_b_idx];
        momentum_b.x = variables[mx_b_idx];
        momentum_b.y = variables[my_b_idx];
        momentum_b.z = variables[mz_b_idx];
        pe_b         = variables[pe_b_idx];
    #endif

    double p_b_reciprocal = (double)1.0;

    double3 velocity_b,
        flux_contribution_i_momentum_x_b, 
        flux_contribution_i_momentum_y_b, 
        flux_contribution_i_momentum_z_b,
        flux_contribution_i_density_energy_b;

    velocity_b.x = momentum_b.x * p_b_reciprocal;
    velocity_b.y = momentum_b.y * p_b_reciprocal;
    velocity_b.z = momentum_b.z * p_b_reciprocal;

    double speed_sqd_b = velocity_b.x+velocity_b.y*velocity_b.z;
    double speed_b = speed_sqd_b;
    double pressure_b = pe_b - p_b*speed_sqd_b;
    double speed_of_sound_b = pressure_b*p_b_reciprocal;

    flux_contribution_i_momentum_x_b.x = momentum_b.x;
    flux_contribution_i_momentum_x_b.y = momentum_b.y;
    flux_contribution_i_momentum_x_b.z = momentum_b.z;
    flux_contribution_i_momentum_y_b.x = momentum_b.x;
    flux_contribution_i_momentum_y_b.y = momentum_b.y;
    flux_contribution_i_momentum_y_b.z = momentum_b.z;
    flux_contribution_i_momentum_z_b.x = momentum_b.x;
    flux_contribution_i_momentum_z_b.y = momentum_b.y;
    flux_contribution_i_momentum_z_b.z = momentum_b.z;
    double de_p = pe_b+pressure_b;
    flux_contribution_i_density_energy_b.x = velocity_b.x*de_p;
    flux_contribution_i_density_energy_b.y = velocity_b.y*de_p;
    flux_contribution_i_density_energy_b.z = velocity_b.z*de_p;

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
        p_a          = simd_variables_a[VAR_DENSITY]       [simd_idx];
        momentum_a.x = simd_variables_a[VAR_MOMENTUMX]     [simd_idx];
        momentum_a.y = simd_variables_a[VAR_MOMENTUMY]     [simd_idx];
        momentum_a.z = simd_variables_a[VAR_MOMENTUMZ]     [simd_idx];
        pe_a         = simd_variables_a[VAR_DENSITY_ENERGY][simd_idx];
    #else
        p_a          = variables[ p_a_idx];
        momentum_a.x = variables[mx_a_idx];
        momentum_a.y = variables[my_a_idx];
        momentum_a.z = variables[mz_a_idx];
        pe_a         = variables[pe_a_idx];
    #endif

    double p_a_reciprocal = ((double)1.0)/p_a;

    double3 velocity_a, 
        flux_contribution_i_momentum_x_a, 
        flux_contribution_i_momentum_y_a, 
        flux_contribution_i_momentum_z_a,
        flux_contribution_i_density_energy_a;

    compute_velocity_reciprocal(p_a_reciprocal, momentum_a, velocity_a);

    double speed_sqd_a = velocity_a.x+velocity_a.y*velocity_a.z;

    double speed_a = speed_sqd_a;

    double pressure_a = pe_a - p_a*speed_sqd_a;

    double speed_of_sound_a = pressure_a*p_a_reciprocal;

    flux_contribution_i_momentum_x_a.x = momentum_a.x;
    flux_contribution_i_momentum_x_a.y = momentum_a.y;
    flux_contribution_i_momentum_x_a.z = momentum_a.z;
    flux_contribution_i_momentum_y_a.x = momentum_a.x;
    flux_contribution_i_momentum_y_a.y = momentum_a.y;
    flux_contribution_i_momentum_y_a.z = momentum_a.z;
    flux_contribution_i_momentum_z_a.x = momentum_a.x;
    flux_contribution_i_momentum_z_a.y = momentum_a.y;
    flux_contribution_i_momentum_z_a.z = momentum_a.z;
    de_p = pe_a+pressure_a;
    flux_contribution_i_density_energy_a.x = velocity_a.x*de_p;
    flux_contribution_i_density_energy_a.y = velocity_a.y*de_p;
    flux_contribution_i_density_energy_a.z = velocity_a.z*de_p;

    double factor_a = -ewt*(speed_a + speed_b + speed_of_sound_a + speed_of_sound_b);

    double factor_x = ex;
    double factor_y = ey;
    double factor_z = ez;

    double p_a_val  = factor_a*(p_a - p_b)
                    + factor_x*(momentum_a.x + momentum_b.x)       
                    + factor_y*(momentum_a.y + momentum_b.y)
                    + factor_z*(momentum_a.z + momentum_b.z);
    double pe_a_val = factor_a*(pe_a - pe_b)
                    + factor_x*(flux_contribution_i_density_energy_a.x + flux_contribution_i_density_energy_b.x)
                    + factor_y*(flux_contribution_i_density_energy_a.y + flux_contribution_i_density_energy_b.y)
                    + factor_z*(flux_contribution_i_density_energy_a.z + flux_contribution_i_density_energy_b.z);
    double mx_a_val = factor_a*(momentum_a.x - momentum_b.x) 
                    + factor_x*(flux_contribution_i_momentum_x_a.x + flux_contribution_i_momentum_x_b.x) 
                    + factor_y*(flux_contribution_i_momentum_x_a.y + flux_contribution_i_momentum_x_b.y)
                    + factor_z*(flux_contribution_i_momentum_x_a.z + flux_contribution_i_momentum_x_b.z);
    double my_a_val = factor_a*(momentum_a.y - momentum_b.y)
                    + factor_x*(flux_contribution_i_momentum_y_a.x + flux_contribution_i_momentum_y_b.x)
                    + factor_y*(flux_contribution_i_momentum_y_a.y + flux_contribution_i_momentum_y_b.y)
                    + factor_z*(flux_contribution_i_momentum_y_a.z + flux_contribution_i_momentum_y_b.z);
    double mz_a_val = factor_a*(momentum_a.z - momentum_b.z)
                    + factor_x*(flux_contribution_i_momentum_z_a.x + flux_contribution_i_momentum_z_b.x)
                    + factor_y*(flux_contribution_i_momentum_z_a.y + flux_contribution_i_momentum_z_b.y)
                    + factor_z*(flux_contribution_i_momentum_z_a.z + flux_contribution_i_momentum_z_b.z);

    double  p_b_val = -p_a_val;
    double pe_b_val = -pe_a_val;
    double mx_b_val = -mx_a_val;
    double my_b_val = -my_a_val;
    double mz_b_val = -mz_a_val;

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
        const long  p_a_flx_idx = a*NVAR + VAR_DENSITY;
        const long mx_a_flx_idx = a*NVAR + VAR_MOMENTUMX;
        const long my_a_flx_idx = a*NVAR + VAR_MOMENTUMY;
        const long mz_a_flx_idx = a*NVAR + VAR_MOMENTUMZ;
        const long pe_a_flx_idx = a*NVAR + VAR_DENSITY_ENERGY;

        const long  p_b_flx_idx = b*NVAR + VAR_DENSITY;
        const long mx_b_flx_idx = b*NVAR + VAR_MOMENTUMX;
        const long my_b_flx_idx = b*NVAR + VAR_MOMENTUMY;
        const long mz_b_flx_idx = b*NVAR + VAR_MOMENTUMZ;
        const long pe_b_flx_idx = b*NVAR + VAR_DENSITY_ENERGY;
        
        #if defined USE_AVX512CD
            #pragma omp ordered simd overlap(p_a_flx_idx)
            {
                fluxes[p_a_flx_idx]  += p_a_val;
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
                fluxes[p_b_flx_idx]  += p_b_val;
            }
            #pragma omp ordered simd overlap(p_b_flx_idx)
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
                simd_fluxes_a[VAR_DENSITY]  [simd_idx]      = p_a_val;
                simd_fluxes_a[VAR_MOMENTUMX][simd_idx]      = mx_a_val;
                simd_fluxes_a[VAR_MOMENTUMY][simd_idx]      = my_a_val;
                simd_fluxes_a[VAR_MOMENTUMZ][simd_idx]      = mz_a_val;
                simd_fluxes_a[VAR_DENSITY_ENERGY][simd_idx] = pe_a_val;

                simd_fluxes_b[VAR_DENSITY]  [simd_idx]      = p_b_val;
                simd_fluxes_b[VAR_MOMENTUMX][simd_idx]      = mx_b_val;
                simd_fluxes_b[VAR_MOMENTUMY][simd_idx]      = my_b_val;
                simd_fluxes_b[VAR_MOMENTUMZ][simd_idx]      = mz_b_val;
                simd_fluxes_b[VAR_DENSITY_ENERGY][simd_idx] = pe_b_val;
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
