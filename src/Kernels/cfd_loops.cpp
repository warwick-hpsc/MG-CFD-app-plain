// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

#include <stdlib.h>

#include "cfd_loops.h"

#include "papi_funcs.h"
#include "timer.h"
#include "loop_stats.h"

// Original step factor calculation from rodinia/cfd
void compute_step_factor_legacy(
    long nel, 
    const double *restrict variables, 
    const double *restrict volumes, 
    double *restrict step_factors)
{
    log("compute_step_factor()");
    current_kernel = COMPUTE_STEP;

    long loop_start = 0;
    long loop_end = nel;

    #ifdef OMP
        #pragma omp parallel firstprivate(loop_start, loop_end)
        {
            openmp_distribute_loop_iterations(&loop_start, &loop_end);
    #endif

    #ifdef PAPI
    start_papi();
    #endif
    #ifdef TIME
    start_timer();
    #endif
    for (long i=loop_start; i<loop_end; i++)
    {
        long p_idx  = NVAR*i + VAR_DENSITY;
        long mx_idx = NVAR*i + VAR_MOMENTUMX;
        long my_idx = NVAR*i + VAR_MOMENTUMY;
        long mz_idx = NVAR*i + VAR_MOMENTUMZ;
        long pe_idx = NVAR*i + VAR_DENSITY_ENERGY;

        double density, density_energy;
        double3 momentum;
        density        = variables[p_idx];
        momentum.x     = variables[mx_idx];
        momentum.y     = variables[my_idx];
        momentum.z     = variables[mz_idx];
        density_energy = variables[pe_idx];

        double3 velocity;
        compute_velocity(density, momentum, velocity);
        double speed_sqd      = compute_speed_sqd(velocity);
        double pressure       = compute_pressure(density, density_energy, speed_sqd);
        double speed_of_sound = compute_speed_of_sound(density, pressure);

        // dt = double(0.5) * std::sqrt(volumes[i]) /  (||v|| + c).... but when we do time stepping, this later would need to be divided by the volume, so we just do it all at once
        step_factors[i] = double(0.5) / (sqrt(volumes[i]) * (sqrt(speed_sqd) + speed_of_sound));
    }
    #ifdef TIME
    stop_timer();
    #endif
    #ifdef PAPI
    stop_papi();
    #endif
    record_iters(loop_start, loop_end);

    #ifdef OMP
    }
    #endif
}

// Corrected step factor calculation
void compute_step_factor(
    long nel, 
    const double *restrict variables, 
    const double *restrict volumes, 
    double *restrict step_factors)
{
    log("compute_step_factor()");
    current_kernel = COMPUTE_STEP;

    long loop_start = 0;
    long loop_end = nel;
    #ifdef OMP
        #pragma omp parallel firstprivate(loop_start, loop_end)
        {
            openmp_distribute_loop_iterations(&loop_start, &loop_end);
    #endif
    #ifdef PAPI
    start_papi();
    #endif
    #ifdef TIME
    start_timer();
    #endif
    for (long i=loop_start; i<loop_end; i++)
    {
        long p_idx  = NVAR*i + VAR_DENSITY;
        long mx_idx = NVAR*i + VAR_MOMENTUMX;
        long my_idx = NVAR*i + VAR_MOMENTUMY;
        long mz_idx = NVAR*i + VAR_MOMENTUMZ;
        long pe_idx = NVAR*i + VAR_DENSITY_ENERGY;

        double density, density_energy;
        double3 momentum;
        density        = variables[p_idx];
        momentum.x     = variables[mx_idx];
        momentum.y     = variables[my_idx];
        momentum.z     = variables[mz_idx];
        density_energy = variables[pe_idx];

        double3 velocity;
        compute_velocity(density, momentum, velocity);
        double speed_sqd      = compute_speed_sqd(velocity);
        double pressure       = compute_pressure(density, density_energy, speed_sqd);
        double speed_of_sound = compute_speed_of_sound(density, pressure);

        // The original step factor calculation in rodinia/cfd refers to 'areas' and  
        // uses square root, implying that calculation is for 2D CFD not 3D. 
        // 'dt' calculation for 3D CFD should use a cubic root, so fix that below:
        double dt = cbrt(volumes[i]) / (sqrt(speed_sqd) + speed_of_sound);
        step_factors[i] = 0.5 * dt;
    }
    #ifdef TIME
    stop_timer();
    #endif
    #ifdef PAPI
    stop_papi();
    #endif
    record_iters(loop_start, loop_end);
    #ifdef OMP
    }
    #endif

    // Sync dt:
    double min_dt = step_factors[0];
    #pragma omp parallel for reduction(min:min_dt)
    for (long i=0; i<nel; i++)
    {
        if (step_factors[i] < min_dt) {
            min_dt = step_factors[i];
        }
    }
    #pragma omp parallel for
    for (long i=0; i<nel; i++)
    {
        step_factors[i] = min_dt;
    }
    // Bring forward a division-by-volume performed by time_step():
    #pragma omp parallel for
    for (long i=0; i<nel; i++)
    {
        step_factors[i] /= volumes[i];
    }
}

void update_edges(
    long first_edge, 
    long nedges, 
    const edge_neighbour *restrict edges,
    const edge *restrict edge_variables,
    double *restrict fluxes)
{
    log("update_edges()");
    current_kernel = UPDATE;

    long loop_start = first_edge;
    long loop_end = loop_start + nedges;

    #if defined OMP && defined OMP_SCATTERS
        #pragma omp parallel firstprivate(loop_start, loop_end)
        {
            openmp_distribute_loop_iterations(&loop_start, &loop_end);
    #endif

    #ifdef PAPI
    start_papi();
    #endif
    #ifdef TIME
    start_timer();
    #endif
    for (long i=loop_start; i<loop_end; i++)
    {
        const long a = edges[i].a;
        if (a >= 0) {
            update_a(
                i, 
                a, fluxes, 
                edge_variables);
        }

        const long b = edges[i].b;
        update_b(
            i,
            b, fluxes, 
            edge_variables);
    }
    #ifdef TIME
    stop_timer();
    #endif
    #ifdef PAPI
    stop_papi();
    #endif
    record_iters(loop_start, loop_end);

    #if defined OMP && defined OMP_SCATTERS
    }
    #endif

    log("update_edges() complete");
}

void time_step(
    int j, 
    long nel, 
    const double *restrict step_factors, 
    double *restrict fluxes, 
    const double *restrict old_variables, 
    double *restrict variables)
{
    log("time_step()");
    current_kernel = TIME_STEP;

    long loop_start = 0;
    long loop_end = nel;

    #ifdef OMP
        #pragma omp parallel firstprivate(loop_start, loop_end)
        {
            openmp_distribute_loop_iterations(&loop_start, &loop_end);
    #endif

    #ifdef PAPI
    start_papi();
    #endif
    #ifdef TIME
    start_timer();
    #endif
    for(long i=loop_start; i<loop_end; i++)
    {
        double factor = (step_factors[i]/double(RK+1-j));

        const long p_var_idx  = NVAR*i + VAR_DENSITY;
        const long mx_var_idx = NVAR*i + VAR_MOMENTUMX;
        const long my_var_idx = NVAR*i + VAR_MOMENTUMY;
        const long mz_var_idx = NVAR*i + VAR_MOMENTUMZ;
        const long pe_var_idx = NVAR*i + VAR_DENSITY_ENERGY;
        
        const long p_flx_idx  = NVAR*i + VAR_DENSITY;
        const long mx_flx_idx = NVAR*i + VAR_MOMENTUMX;
        const long my_flx_idx = NVAR*i + VAR_MOMENTUMY;
        const long mz_flx_idx = NVAR*i + VAR_MOMENTUMZ;
        const long pe_flx_idx = NVAR*i + VAR_DENSITY_ENERGY;

        variables[p_var_idx]  = old_variables[p_var_idx]  + factor*fluxes[p_flx_idx];
        variables[mx_var_idx] = old_variables[mx_var_idx] + factor*fluxes[mx_flx_idx];
        variables[my_var_idx] = old_variables[my_var_idx] + factor*fluxes[my_flx_idx];
        variables[mz_var_idx] = old_variables[mz_var_idx] + factor*fluxes[mz_flx_idx];
        variables[pe_var_idx] = old_variables[pe_var_idx] + factor*fluxes[pe_flx_idx];

        fluxes[p_flx_idx]  = 0.0;
        fluxes[mx_flx_idx] = 0.0;
        fluxes[my_flx_idx] = 0.0;
        fluxes[mz_flx_idx] = 0.0;
        fluxes[pe_flx_idx] = 0.0;
    }
    #ifdef TIME
    stop_timer();
    #endif
    #ifdef PAPI
    stop_papi();
    #endif
    record_iters(loop_start, loop_end);

    #ifdef OMP
    }
    #endif
}

void zero_fluxes(
    long nel, 
    double *restrict fluxes)
{
    log("zero_fluxes()");

    #ifdef OMP
        #pragma omp parallel for
    #endif
    for(long i=0; i<nel; i++)
    {
        const long p_idx  = NVAR*i + VAR_DENSITY;
        const long mx_idx = NVAR*i + VAR_MOMENTUMX;
        const long my_idx = NVAR*i + VAR_MOMENTUMY;
        const long mz_idx = NVAR*i + VAR_MOMENTUMZ;
        const long pe_idx = NVAR*i + VAR_DENSITY_ENERGY;

        fluxes[p_idx]  = 0.0;
        fluxes[mx_idx] = 0.0;
        fluxes[my_idx] = 0.0;
        fluxes[mz_idx] = 0.0;
        fluxes[pe_idx] = 0.0;
    }
}
