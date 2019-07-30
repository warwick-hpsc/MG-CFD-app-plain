#include "flux_loops.h"
#include "cfd_loops.h"

#include "flux_kernel.h"
#ifdef FLUX_CRIPPLE
    #include "flux_kernel_crippled.h"
#endif

#include "papi_funcs.h"
#include "timer.h"
#include "loop_stats.h"

#include <omp.h>

void compute_boundary_flux_edge(
    int first_edge,
    int nedges, 
    const edge_neighbour *restrict edges, 
    const double *restrict variables, 
    #ifdef FLUX_FISSION
        edge *restrict edge_variables
    #else
        double *restrict fluxes
    #endif
    )
{
    log("Computing boundary flux");
    current_kernel = COMPUTE_FLUX_EDGE;

    int loop_start = first_edge;
    int loop_end = loop_start + nedges;

    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        #pragma omp parallel firstprivate(loop_start, loop_end)
        {
            openmp_distribute_loop_iterations(&loop_start, &loop_end);
    #endif
    for (int i=loop_start; i<loop_end; i++)
    {
        #include "flux_boundary_kernel.elemfunc.c"
    }
    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        }
    #endif

    log("Boundary flux compute complete");
}

void compute_wall_flux_edge(
    int first_edge,
    int nedges,
    const edge_neighbour *restrict edges, 
    const double *restrict variables, 
    #ifdef FLUX_FISSION
        edge *restrict edge_variables
    #else
        double *restrict fluxes
    #endif
    )
{
    log("Computing wall flux");
    current_kernel = COMPUTE_FLUX_EDGE;

    int loop_start = first_edge;
    int loop_end = loop_start + nedges;

    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        #pragma omp parallel firstprivate(loop_start, loop_end)
        {
            openmp_distribute_loop_iterations(&loop_start, &loop_end);
    #endif
    for (int i=loop_start; i<loop_end; i++)
    {
        #include "flux_wall_kernel.elemfunc.c"
    }
    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        }
    #endif

    log("Wall flux compute complete");
}

void compute_flux_edge(
    int first_edge,
    int nedges,
    const edge_neighbour *restrict edges, 
    #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
        const double *restrict edge_weights,
    #endif
    const double *restrict variables, 
    #ifdef FLUX_FISSION
        edge *restrict edge_variables
    #else
        double *restrict fluxes
        #ifdef COLOURED_CONFLICT_AVOIDANCE
        , int serial_section_start
        #endif
    #endif
    )
{
    log("Computing internal flux");
    current_kernel = COMPUTE_FLUX_EDGE;

    int loop_start = first_edge;
    int loop_end = loop_start + nedges;
    #if defined SIMD && defined COLOURED_CONFLICT_AVOIDANCE
        loop_end = serial_section_start;
    #endif

    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        #pragma omp parallel firstprivate(loop_start, loop_end)
        {
            openmp_distribute_loop_iterations(&loop_start, &loop_end);
    #endif

    #ifdef FLUX_CRIPPLE
        iters_monitoring_state = 0;
    #else
        #ifdef PAPI
        start_papi();
        #endif
        #ifdef TIME
        start_timer();
        #endif
    #endif

    #ifndef SIMD
        #pragma omp simd safelen(1)
    #else
        #ifdef FLUX_FISSION
            // SIMD is safe
            #pragma omp simd simdlen(DBLS_PER_SIMD)
        #else
            // Conflict avoidance is required for safe SIMD
            #if defined __AVX512CD__ && defined __ICC
                #pragma omp simd safelen(1)
                // TODO: Insert the following pragma into flux_kernel.elemfunc.c to 
                //       enable safe AVX-512-CD SIMD:
                // #pragma omp ordered simd overlap(...)
            #elif defined COLOURED_CONFLICT_AVOIDANCE
                #pragma omp simd simdlen(DBLS_PER_SIMD)
            #endif
        #endif
    #endif
    for (int i=loop_start; i<loop_end; i++)
    {
        compute_flux_edge_kernel(
            #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                edge_weights[i],
            #endif
            edges[i].x, 
            edges[i].y, 
            edges[i].z,
            &variables[edges[i].a*NVAR],
            &variables[edges[i].b*NVAR],
            #ifdef FLUX_FISSION
                &edge_variables[i]
            #else
                &fluxes[edges[i].a*NVAR],
                &fluxes[edges[i].b*NVAR]
            #endif
            );
    }
    #ifdef FLUX_CRIPPLE
        iters_monitoring_state = 1;
    #else
        #ifdef TIME
        stop_timer();
        #endif
        #ifdef PAPI
        stop_papi();
        #endif
        record_iters(loop_start, loop_end);
    #endif

    #if defined SIMD && defined COLOURED_CONFLICT_AVOIDANCE
        // Compute fluxes of 'remainder' edges without SIMD:
        #ifdef COLOURED_CONFLICT_AVOIDANCE
            loop_start = serial_section_start;
            loop_end = first_edge + nedges;
            #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
                // All remainder edges are separate from the SIMD-able edges, so 
                // a workload distribution is required.
                openmp_distribute_loop_iterations(&loop_start, &loop_end);
            #endif
        #endif

        #pragma omp simd safelen(1)
        for (int i=loop_start; i<loop_end; i++)
        {
            compute_flux_edge_kernel(
                #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                    edge_weights[i],
                #endif
                edges[i].x, 
                edges[i].y, 
                edges[i].z,
                &variables[edges[i].a*NVAR],
                &variables[edges[i].b*NVAR],
                #ifdef FLUX_FISSION
                    &edge_variables[i]
                #else
                    #ifdef MANUAL_CONFLICT_AVOIDANCE
                        i-loop_start,
                        fluxes_a, 
                        fluxes_b
                    #else
                        &fluxes[edges[i].a*NVAR],
                        &fluxes[edges[i].b*NVAR]
                    #endif
                #endif
                );
        }
    #endif

    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        }
    #endif

    log("Internal flux compute complete");
}

#ifdef FLUX_CRIPPLE
void compute_flux_edge_crippled(
    int first_edge,
    int nedges,
    const edge_neighbour *restrict edges, 
    #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
        const double *restrict edge_weights,
    #endif
    const double *restrict variables, 
    #ifdef FLUX_FISSION
        edge *restrict edge_variables
    #else
        double *restrict fluxes
        #ifdef COLOURED_CONFLICT_AVOIDANCE
        , int serial_section_start
        #endif
    #endif
    )
{
    log("Computing internal flux crippled");
    current_kernel = COMPUTE_FLUX_EDGE;

    int loop_start = first_edge;
    int loop_end = loop_start + nedges;
    #if defined SIMD && defined COLOURED_CONFLICT_AVOIDANCE
        loop_end = serial_section_start;
    #endif

    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
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

    #ifndef SIMD
        #pragma omp simd safelen(1)
    #else
        #ifdef FLUX_FISSION
            // SIMD is safe
            #pragma omp simd simdlen(DBLS_PER_SIMD)
        #else
            // Conflict avoidance is required for safe SIMD
            #if defined __AVX512CD__ && defined __ICC
                #pragma omp simd safelen(1)
                // TODO: Insert the following pragma into flux.elemfunc.c to 
                //       enable safe AVX-512-CD SIMD:
                // #pragma omp ordered simd overlap(...)
            #elif defined COLOURED_CONFLICT_AVOIDANCE
                #pragma omp simd simdlen(DBLS_PER_SIMD)
            #endif
        #endif
    #endif
    for (int i=loop_start; i<loop_end; i++)
    {
        compute_flux_edge_kernel_crippled(
            #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                edge_weights[i],
            #else
                edges[i].x, 
                edges[i].y, 
                edges[i].z,
            #endif
            &variables[edges[i].a*NVAR],
            &variables[edges[i].b*NVAR],
            #ifdef FLUX_FISSION
                &edge_variables[i]
            #else
                &fluxes[edges[i].a*NVAR],
                &fluxes[edges[i].b*NVAR]
            #endif
            );
    }

    #if defined SIMD && COLOURED_CONFLICT_AVOIDANCE
        // Compute fluxes of 'remainder' edges without SIMD:
        #ifdef COLOURED_CONFLICT_AVOIDANCE
            loop_start = serial_section_start;
            loop_end = first_edge + nedges;
            #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
                // All remainder edges are separate from the SIMD-able edges, so 
                // a workload distribution is required.
                openmp_distribute_loop_iterations(&loop_start, &loop_end);
            #endif
        #endif
        
        #pragma omp simd safelen(1)
        for (int i=loop_start; i<loop_end; i++)
        {
            compute_flux_edge_kernel(
                #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                    edge_weights[i],
                #else
                    edges[i].x, 
                    edges[i].y, 
                    edges[i].z,
                #endif
                &variables[edges[i].a*NVAR],
                &variables[edges[i].b*NVAR],
                #ifdef FLUX_FISSION
                    &edge_variables[i]
                #else
                    #ifdef MANUAL_CONFLICT_AVOIDANCE
                        i-loop_start,
                        fluxes_a, 
                        fluxes_b
                    #else
                        &fluxes[edges[i].a*NVAR],
                        &fluxes[edges[i].b*NVAR]
                    #endif
                #endif
                );
        }
    #endif

    #ifdef TIME
    stop_timer();
    #endif
    #ifdef PAPI
    stop_papi();
    #endif
    record_iters(loop_start, loop_end);

    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        }
    #endif
}
#endif
