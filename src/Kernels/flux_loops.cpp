#include "flux_loops.h"
#include "cfd_loops.h"

#include "flux_kernel.h"
#ifdef FLUX_CRIPPLE
    #include "flux_kernel_crippled.h"
#endif

#include "papi_funcs.h"
#include "timer.h"
#include "loop_stats.h"

void compute_boundary_flux_edge(
    long first_edge,
    long nedges, 
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

    long loop_start = first_edge;
    long loop_end = loop_start + nedges;

    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        #pragma omp parallel firstprivate(loop_start, loop_end)
        {
            openmp_distribute_loop_iterations(&loop_start, &loop_end);
    #endif
    #ifndef SIMD
        #pragma omp simd safelen(1)
    #endif
    for (long i=loop_start; i<loop_end; i++)
    {
        #include "flux_boundary_kernel.elemfunc.c"
    }
    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        }
    #endif

    log("Boundary flux compute complete");
}

void compute_wall_flux_edge(
    long first_edge,
    long nedges,
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

    long loop_start = first_edge;
    long loop_end = loop_start + nedges;

    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        #pragma omp parallel firstprivate(loop_start, loop_end)
        {
            openmp_distribute_loop_iterations(&loop_start, &loop_end);
    #endif
    #ifndef SIMD
        #pragma omp simd safelen(1)
    #endif
    for (long i=loop_start; i<loop_end; i++)
    {
        #include "flux_wall_kernel.elemfunc.c"
    }
    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        }
    #endif

    log("Wall flux compute complete");
}

void compute_flux_edge(
    long first_edge,
    long nedges,
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

    long loop_start = first_edge;
    long loop_end = loop_start + nedges;
    #if defined SIMD && defined COLOURED_CONFLICT_AVOIDANCE && (!defined FLUX_FISSION)
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
        record_iters(loop_start, loop_end);
    #endif

    #ifndef SIMD
        #pragma omp simd safelen(1)
    #else
        #ifdef FLUX_FISSION
            // SIMD is safe
            #pragma omp simd simdlen(DBLS_PER_SIMD)
        #else
            // Conflict avoidance is required for safe SIMD
            #if defined COLOURED_CONFLICT_AVOIDANCE
                #pragma omp simd simdlen(DBLS_PER_SIMD)
            #elif defined MANUAL_CONFLICT_AVOIDANCE
                const int loop_start_orig = loop_start;
                const int loop_end_orig = loop_end;
                int v_start = loop_start;
                int v_end = loop_start + ((loop_end-loop_start)/DBLS_PER_SIMD)*DBLS_PER_SIMD;
                double fluxes_a[NVAR][DBLS_PER_SIMD];
                double fluxes_b[NVAR][DBLS_PER_SIMD];
                for (int v=v_start; v<v_end; v+=DBLS_PER_SIMD) {
                    for (int x=0; x<NVAR; x++) {
                        for (int n=0; n<DBLS_PER_SIMD; n++) {
                            fluxes_a[x][n] = 0.0;
                            fluxes_b[x][n] = 0.0;
                        }
                    }
                    loop_start = v;
                    loop_end = v+DBLS_PER_SIMD;
                    #pragma omp simd simdlen(DBLS_PER_SIMD)

            #elif defined USE_AVX512CD
                // Always prefer using OMP pragma to vectorise, gives better performance 
                // than default auto-vectoriser triggered by absent pragma
                #pragma omp simd simdlen(DBLS_PER_SIMD)
            #endif
        #endif
    #endif
    for (long i=loop_start; i<loop_end; i++)
    {
        #if defined USE_AVX512CD
            // For Intel AVX-512-CD auto-vectorizer to act, I need to 
            // directly include the kernel source here rather than 
            // call an inlined kernel function.
            #include "flux_kernel.elemfunc.c"
        #else
            compute_flux_edge_kernel(
                edges[i].a, edges[i].b, 
                #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                    edge_weights[i],
                #endif
                edges[i].x, 
                edges[i].y, 
                edges[i].z,
                &variables[edges[i].a*NVAR],
                &variables[edges[i].b*NVAR],
                #ifdef FLUX_FISSION
                    &edge_variables[i*NVAR]
                #else
                    #if defined SIMD and defined MANUAL_CONFLICT_AVOIDANCE
                        i-loop_start,
                        fluxes_a, 
                        fluxes_b
                    #else
                        &fluxes[edges[i].a*NVAR],
                        &fluxes[edges[i].b*NVAR]
                    #endif
                #endif
                );
        #endif
    }

    #if defined SIMD && defined MANUAL_CONFLICT_AVOIDANCE && (!defined FLUX_FISSION)
        // Write out fluxes:
            for (int x=0; x<NVAR; x++) {
                #pragma omp simd safelen(1)
                for (int n=0; n<DBLS_PER_SIMD; n++) {
                    int a = edges[v+n].a;
                    int b = edges[v+n].b;
                    fluxes[a*NVAR+x] += fluxes_a[x][n];
                    fluxes[b*NVAR+x] += fluxes_b[x][n];
                }
            }
        } // Close outer loop over SIMD blocks
    #endif

    #ifdef FLUX_CRIPPLE
        iters_monitoring_state = 1;
    #else
        #ifdef TIME
        stop_timer();
        #endif
        #ifdef PAPI
        stop_papi();
        #endif
    #endif

    #if defined SIMD && (defined COLOURED_CONFLICT_AVOIDANCE || defined MANUAL_CONFLICT_AVOIDANCE) && (!defined FLUX_FISSION)
        // Compute fluxes of 'remainder' edges without SIMD:
        #ifdef COLOURED_CONFLICT_AVOIDANCE
            loop_start = serial_section_start;
            loop_end = first_edge + nedges;
            #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
                // All remainder edges are separate from the SIMD-able edges, so 
                // a workload distribution is required.
                openmp_distribute_loop_iterations(&loop_start, &loop_end);
            #endif
        #elif defined MANUAL_CONFLICT_AVOIDANCE
            for (int x=0; x<NVAR; x++) {
                for (int i=0; i<DBLS_PER_SIMD; i++) {
                    fluxes_a[x][i] = 0.0;
                    fluxes_b[x][i] = 0.0;
                }
            }
            loop_start = v_end;
            loop_end = loop_end_orig;
            #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
                // Each thread already has its own remainders to process, no need to 
                // distribute workload
            #endif
        #endif

        #pragma omp simd safelen(1)
        for (int i=loop_start; i<loop_end; i++)
        {
            compute_flux_edge_kernel(
                edges[i].a, edges[i].b,
                #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                    edge_weights[i],
                #endif
                edges[i].x, 
                edges[i].y, 
                edges[i].z,
                &variables[edges[i].a*NVAR],
                &variables[edges[i].b*NVAR],
                #ifdef FLUX_FISSION
                    &edge_variables[i*NVAR]
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
        #ifdef MANUAL_CONFLICT_AVOIDANCE
            // Write out fluxes:
            #pragma omp simd safelen(1)
            for (int i=loop_start; i<loop_end; i++) {
                int a = edges[i].a;
                int b = edges[i].b;
                for (int x=0; x<NVAR; x++) {
                    fluxes[a*NVAR+x] += fluxes_a[x][i-loop_start];
                    fluxes[b*NVAR+x] += fluxes_b[x][i-loop_start];
                }
            }
        #endif
    #endif

    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        }
    #endif

    log("Internal flux compute complete");
}

#ifdef FLUX_CRIPPLE
void compute_flux_edge_crippled(
    long first_edge,
    long nedges,
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

    long loop_start = first_edge;
    long loop_end = loop_start + nedges;
    #if defined SIMD && defined COLOURED_CONFLICT_AVOIDANCE && (!defined FLUX_FISSION)
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
    record_iters(loop_start, loop_end);

    #ifndef SIMD
        #pragma omp simd safelen(1)
    #else
        #ifdef FLUX_FISSION
            // SIMD is safe
            #pragma omp simd simdlen(DBLS_PER_SIMD)
        #else
            // Conflict avoidance is required for safe SIMD
            #if defined COLOURED_CONFLICT_AVOIDANCE
                #pragma omp simd simdlen(DBLS_PER_SIMD)
            #elif defined MANUAL_CONFLICT_AVOIDANCE
                const int loop_start_orig = loop_start;
                const int loop_end_orig = loop_end;
                int v_start = loop_start;
                int v_end = loop_start + ((loop_end-loop_start)/DBLS_PER_SIMD)*DBLS_PER_SIMD;
                double fluxes_a[NVAR][DBLS_PER_SIMD];
                double fluxes_b[NVAR][DBLS_PER_SIMD];
                for (int v=v_start; v<v_end; v+=DBLS_PER_SIMD) {
                    for (int x=0; x<NVAR; x++) {
                        for (int n=0; n<DBLS_PER_SIMD; n++) {
                            fluxes_a[x][n] = 0.0;
                            fluxes_b[x][n] = 0.0;
                        }
                    }
                    loop_start = v;
                    loop_end = v+DBLS_PER_SIMD;
                    #pragma omp simd simdlen(DBLS_PER_SIMD)

            #elif defined USE_AVX512CD
                // Always prefer using OMP pragma to vectorise, gives better performance 
                // than default auto-vectoriser triggered by absent pragma
                #pragma omp simd simdlen(DBLS_PER_SIMD)
            #endif
        #endif
    #endif
    for (long i=loop_start; i<loop_end; i++)
    {
        #if defined USE_AVX512CD
            // For Intel AVX-512-CD auto-vectorizer to act, I need to 
            // directly include the kernel source here rather than 
            // call an inlined kernel function.
            #include "flux_kernel_crippled.elemfunc.c"
        #else
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
                    &edge_variables[i*NVAR]
                #else
                    #if defined SIMD and defined MANUAL_CONFLICT_AVOIDANCE
                        i-loop_start,
                        fluxes_a, 
                        fluxes_b
                    #else
                        &fluxes[edges[i].a*NVAR],
                        &fluxes[edges[i].b*NVAR]
                    #endif
                #endif
                );
        #endif
    }

    #if defined SIMD && defined MANUAL_CONFLICT_AVOIDANCE && (!defined FLUX_FISSION)
        // Write out fluxes:
            for (int x=0; x<NVAR; x++) {
                #pragma omp simd safelen(1)
                for (int n=0; n<DBLS_PER_SIMD; n++) {
                    int a = edges[v+n].a;
                    int b = edges[v+n].b;
                    fluxes[a*NVAR+x] += fluxes_a[x][n];
                    fluxes[b*NVAR+x] += fluxes_b[x][n];
                }
            }
        } // Close outer loop over SIMD blocks
    #endif

    #if defined SIMD && (defined COLOURED_CONFLICT_AVOIDANCE || defined MANUAL_CONFLICT_AVOIDANCE) && (!defined FLUX_FISSION)
        // Compute fluxes of 'remainder' edges without SIMD:
        #ifdef COLOURED_CONFLICT_AVOIDANCE
            loop_start = serial_section_start;
            loop_end = first_edge + nedges;
            #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
                // All remainder edges are separate from the SIMD-able edges, so 
                // a workload distribution is required.
                openmp_distribute_loop_iterations(&loop_start, &loop_end);
            #endif
        #elif defined MANUAL_CONFLICT_AVOIDANCE
            for (int i=0; i<DBLS_PER_SIMD; i++) {
                for (int x=0; x<NVAR; x++) {
                    fluxes_a[x][i] = 0.0;
                    fluxes_b[x][i] = 0.0;
                }
            }
            loop_start = v_end;
            loop_end = loop_end_orig;
            #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
                // Each thread already has its own remainders to process, no need to 
                // distribute workload
            #endif
        #endif
        
        #pragma omp simd safelen(1)
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
                    &edge_variables[i*NVAR]
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
        #ifdef MANUAL_CONFLICT_AVOIDANCE
            // Write out fluxes:
            #pragma omp simd safelen(1)
            for (int i=loop_start; i<loop_end; i++) {
                for (int v=0; v<NVAR; v++) {
                    fluxes[edges[i].a*NVAR+v] += fluxes_a[v][i-loop_start];
                    fluxes[edges[i].b*NVAR+v] += fluxes_b[v][i-loop_start];
                }
            }
        #endif
    #endif

    #ifdef TIME
    stop_timer();
    #endif
    #ifdef PAPI
    stop_papi();
    #endif

    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        }
    #endif
}
#endif
