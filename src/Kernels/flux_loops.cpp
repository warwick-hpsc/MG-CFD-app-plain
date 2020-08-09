#include "flux_loops.h"
#include "cfd_loops.h"

#include "flux_kernel.h"
#ifdef FLUX_CRIPPLE
    #include "flux_kernel_crippled.h"
#endif

#include "papi_funcs.h"
#include "timer.h"
#include "loop_stats.h"

void compute_flux_edge(
    long first_edge,
    long nedges,
    const edge_neighbour *restrict edges, 
    // edge_neighbour *restrict edges, 
    #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
        const double *restrict edge_weights,
    #endif
    const double *restrict variables, 
    #ifdef FLUX_FISSION
        edge *restrict edge_variables
    #else
        double *restrict fluxes
        #ifdef COLOURED_CONFLICT_AVOIDANCE
        , long serial_section_start
        #endif
    #endif
    )
{
    log("Computing internal flux");
    current_kernel = COMPUTE_FLUX_EDGE;

    long flux_loop_start = first_edge;
    long loop_end = flux_loop_start + nedges;

    // edges = (edge_neighbour*)__builtin_assume_aligned(edges, 64);

    #if defined SIMD && defined COLOURED_CONFLICT_AVOIDANCE && (!defined FLUX_FISSION)
        loop_end = serial_section_start;
    #endif

    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        #pragma omp parallel firstprivate(flux_loop_start, loop_end)
        {
            openmp_distribute_loop_iterations(&flux_loop_start, &loop_end);
    #endif

    #if defined MANUAL_SCATTER
        double fluxes_a[NVAR][MANUAL_WIDTH];
        double fluxes_b[NVAR][MANUAL_WIDTH];
        #ifdef MANUAL_GATHER
            double variables_a[NVAR][MANUAL_WIDTH];
            double variables_b[NVAR][MANUAL_WIDTH];
        #endif
        // double fluxes_a[NVAR][MANUAL_WIDTH] __attribute__((aligned(64)));
        // double fluxes_b[NVAR][MANUAL_WIDTH] __attribute__((aligned(64)));
        // #ifdef MANUAL_GATHER
        //     double variables_a[NVAR][MANUAL_WIDTH] __attribute__((aligned(64)));
        //     double variables_b[NVAR][MANUAL_WIDTH] __attribute__((aligned(64)));
        // #endif
        for (int n=0; n<MANUAL_WIDTH; n++) {
            for (int x=0; x<NVAR; x++) {
                fluxes_a[x][n] = 0.0;
                fluxes_b[x][n] = 0.0;
            }
        }
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
        record_iters(flux_loop_start, loop_end);
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
                // Preventing unrolling of outer loop helps assembly-loop-extractor
                #pragma nounroll
            #elif defined MANUAL_SCATTER
                const long loop_end_orig = loop_end;
                long v_start = flux_loop_start;
                long v_end = flux_loop_start + ((loop_end-flux_loop_start)/MANUAL_WIDTH)*MANUAL_WIDTH;
                for (long v=v_start; v<v_end; v+=MANUAL_WIDTH) {
                    #ifdef MANUAL_GATHER
                        #pragma omp simd safelen(1)
                        // Preventing unrolling of outer loop helps assembly-loop-extractor
                        #pragma nounroll
                        for (int n=0; n<MANUAL_WIDTH; n++) {
                            #pragma unroll
                            for (int x=0; x<NVAR; x++) {
                                variables_a[x][n] = variables[edges[v+n].a*NVAR+x];
                                variables_b[x][n] = variables[edges[v+n].b*NVAR+x];
                            }
                        }
                    #endif
                    flux_loop_start = v;
                    loop_end = v+MANUAL_WIDTH;

                    #ifdef __clang__
                        // __builtin_assume(flux_loop_start%DBLS_PER_SIMD == 0);
                        // __builtin_assume(loop_end%DBLS_PER_SIMD == 0);
                    #endif
                    #pragma omp simd simdlen(DBLS_PER_SIMD)

            #elif defined USE_AVX512CD
                // Always prefer using OMP pragma to vectorise, gives better performance 
                // than default auto-vectoriser triggered by absent pragma
                #pragma omp simd simdlen(DBLS_PER_SIMD)
            #endif
        #endif
    #endif
    for (long i=flux_loop_start; i<loop_end; i++)
    {
        #if defined USE_AVX512CD
            // For Intel AVX-512-CD auto-vectorizer to act, I need to 
            // directly include the kernel source here rather than 
            // call an inlined kernel function.
            #include "flux_kernel.elemfunc.c"
        #elif defined __clang__
            // Clang is not inlining call to compute_flux_edge_kernel(), 
            // so manually inline here:
            #include "flux_kernel.elemfunc.c"
        #elif defined __GNUC__
            // Inlining call makes GNU vector log easier to parse
            #include "flux_kernel.elemfunc.c"
        #else
            compute_flux_edge_kernel(
                #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                    edge_weights[i],
                #endif
                edges[i].x, 
                edges[i].y, 
                edges[i].z,
                #if defined SIMD && defined MANUAL_GATHER
                    variables_a, 
                    variables_b,
                #else
                    &variables[edges[i].a*NVAR],
                    &variables[edges[i].b*NVAR],
                #endif
                #ifdef FLUX_FISSION
                    &edge_variables[i*NVAR]
                #else
                    #if defined SIMD && defined MANUAL_SCATTER
                        i-flux_loop_start,
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

    #if defined SIMD && defined MANUAL_SCATTER && (!defined FLUX_FISSION)
        // Write out fluxes:
            // Preventing unrolling of outer loop helps assembly-loop-extractor
            #pragma nounroll
            for (long n=0; n<MANUAL_WIDTH; n++) {
                #pragma omp simd safelen(1)
                #pragma unroll
                for (long x=0; x<NVAR; x++) {
                    fluxes[edges[v+n].a*NVAR+x] += fluxes_a[x][n];
                    fluxes[edges[v+n].b*NVAR+x] += fluxes_b[x][n];
                    fluxes_a[x][n] = 0.0;
                    fluxes_b[x][n] = 0.0;
                    // Note: zeroing fluxes[] here prevents Clang replacing 
                    // with LLVM memset() calls, which prevents assembly-loop-extractor 
                    // confidently identifying compute loop.
                }
            }
        } // Close outer loop over SIMD blocks
    #endif

    #if defined SIMD && (defined COLOURED_CONFLICT_AVOIDANCE || defined MANUAL_SCATTER) && (!defined FLUX_FISSION)
        // Compute fluxes of 'remainder' edges without SIMD:
        #ifdef COLOURED_CONFLICT_AVOIDANCE
            long remainder_loop_start = serial_section_start;
            loop_end = first_edge + nedges;
            #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
                // All remainder edges are separate from the SIMD-able edges, so 
                // a workload distribution is required.
                openmp_distribute_loop_iterations(&remainder_loop_start, &loop_end);
            #endif
        #elif defined MANUAL_SCATTER
            long remainder_loop_start = v_end;
            loop_end = loop_end_orig;
            for (int x=0; x<NVAR; x++) {
                for (long i=remainder_loop_start; i<loop_end; i++) {
                    fluxes_a[x][i-remainder_loop_start] = 0.0;
                    fluxes_b[x][i-remainder_loop_start] = 0.0;
                    #ifdef MANUAL_GATHER
                        variables_a[x][i-remainder_loop_start] = variables[edges[i].a*NVAR+x];
                        variables_b[x][i-remainder_loop_start] = variables[edges[i].b*NVAR+x];
                    #endif
                }
            }
            #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
                // Each thread already has its own remainders to process, no need to 
                // distribute workload
            #endif
        #endif

        #pragma omp simd safelen(1)
        for (long i=remainder_loop_start; i<loop_end; i++)
        {
            compute_flux_edge_kernel(
                #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                    edge_weights[i],
                #endif
                edges[i].x, 
                edges[i].y, 
                edges[i].z,
                #if defined SIMD && defined MANUAL_GATHER
                    variables_a, 
                    variables_b,
                #else
                    &variables[edges[i].a*NVAR],
                    &variables[edges[i].b*NVAR],
                #endif
                #ifdef FLUX_FISSION
                    &edge_variables[i*NVAR]
                #else
                    #ifdef MANUAL_SCATTER
                        i-remainder_loop_start,
                        fluxes_a, 
                        fluxes_b
                    #else
                        &fluxes[edges[i].a*NVAR],
                        &fluxes[edges[i].b*NVAR]
                    #endif
                #endif
                );
        }
        #ifdef MANUAL_SCATTER
            // Write out fluxes:
            for (long i=remainder_loop_start; i<loop_end; i++) {
                for (int x=0; x<NVAR; x++) {
                    fluxes[edges[i].a*NVAR+x] += fluxes_a[x][i-remainder_loop_start];
                    fluxes[edges[i].b*NVAR+x] += fluxes_b[x][i-remainder_loop_start];
                }
            }
        #endif
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
        , long serial_section_start
        #endif
    #endif
    )
{
    log("Computing internal flux crippled");
    current_kernel = COMPUTE_FLUX_EDGE;

    long flux_loop_start = first_edge;
    long loop_end = flux_loop_start + nedges;
    #if defined SIMD && defined COLOURED_CONFLICT_AVOIDANCE && (!defined FLUX_FISSION)
        loop_end = serial_section_start;
    #endif

    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        #pragma omp parallel firstprivate(flux_loop_start, loop_end)
        {
            openmp_distribute_loop_iterations(&flux_loop_start, &loop_end);
    #endif

    #if defined MANUAL_SCATTER
        double fluxes_a[NVAR][MANUAL_WIDTH];
        double fluxes_b[NVAR][MANUAL_WIDTH];
        #ifdef MANUAL_GATHER
            double variables_a[NVAR][MANUAL_WIDTH];
            double variables_b[NVAR][MANUAL_WIDTH];
        #endif
        for (int n=0; n<MANUAL_WIDTH; n++) {
            for (int x=0; x<NVAR; x++) {
                fluxes_a[x][n] = 0.0;
                fluxes_b[x][n] = 0.0;
            }
        }
    #endif

    #ifdef PAPI
    start_papi();
    #endif
    #ifdef TIME
    start_timer();
    #endif
    record_iters(flux_loop_start, loop_end);

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
                #pragma nounroll
            #elif defined MANUAL_SCATTER
                const long loop_end_orig = loop_end;
                long v_start = flux_loop_start;
                long v_end = flux_loop_start + ((loop_end-flux_loop_start)/MANUAL_WIDTH)*MANUAL_WIDTH;
                #pragma nounroll
                for (long v=v_start; v<v_end; v+=MANUAL_WIDTH) {
                    #ifdef MANUAL_GATHER
                        // Preventing unrolling of outer loop helps assembly-loop-extractor
                        #pragma nounroll
                        #pragma omp simd safelen(1)
                        for (int n=0; n<MANUAL_WIDTH; n++) {
                            #pragma unroll
                            for (int x=0; x<NVAR; x++) {
                                variables_a[x][n] = variables[edges[v+n].a*NVAR+x];
                                variables_b[x][n] = variables[edges[v+n].b*NVAR+x];
                            }
                        }
                    #endif
                    flux_loop_start = v;
                    loop_end = v+MANUAL_WIDTH;

                    #pragma omp simd simdlen(DBLS_PER_SIMD)

            #elif defined USE_AVX512CD
                // Always prefer using OMP pragma to vectorise, gives better performance 
                // than default auto-vectoriser triggered by absent pragma
                #pragma omp simd simdlen(DBLS_PER_SIMD)
            #endif
        #endif
    #endif
    for (long i=flux_loop_start; i<loop_end; i++)
    {
        #if defined USE_AVX512CD
            // For Intel AVX-512-CD auto-vectorizer to act, I need to 
            // directly include the kernel source here rather than 
            // call an inlined kernel function.
            #include "flux_kernel_crippled.elemfunc.c"
        #elif defined __clang__
            // Clang is not inlining call, so manually inline here:
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
                #if defined SIMD && defined MANUAL_GATHER
                    variables_a, 
                    variables_b,
                #else
                    &variables[edges[i].a*NVAR],
                    &variables[edges[i].b*NVAR],
                #endif
                #ifdef FLUX_FISSION
                    &edge_variables[i*NVAR]
                #else
                    #if defined SIMD && defined MANUAL_SCATTER
                        i-flux_loop_start,
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

    #if defined SIMD && defined MANUAL_SCATTER && (!defined FLUX_FISSION)
        // Write out fluxes:
            // Preventing unrolling of outer loop helps assembly-loop-extractor
            #pragma nounroll
            for (long n=0; n<MANUAL_WIDTH; n++) {
                #pragma unroll
                #pragma omp simd safelen(1)
                for (long x=0; x<NVAR; x++) {
                    fluxes[edges[v+n].a*NVAR+x] += fluxes_a[x][n];
                    fluxes[edges[v+n].b*NVAR+x] += fluxes_b[x][n];
                    fluxes_a[x][n] = 0.0;
                    fluxes_b[x][n] = 0.0;
                    // Note: zeroing fluxes[] here prevents Clang replacing 
                    // with LLVM memset() calls, which prevents assembly-loop-extractor 
                    // confidently identifying compute loop.
                }
            }
        } // Close outer loop over SIMD blocks
    #endif

    #if defined SIMD && (defined COLOURED_CONFLICT_AVOIDANCE || defined MANUAL_SCATTER) && (!defined FLUX_FISSION)
        // Compute fluxes of 'remainder' edges without SIMD:
        #ifdef COLOURED_CONFLICT_AVOIDANCE
            long remainder_loop_start = serial_section_start;
            loop_end = first_edge + nedges;
            #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
                // All remainder edges are separate from the SIMD-able edges, so 
                // a workload distribution is required.
                openmp_distribute_loop_iterations(&remainder_loop_start, &loop_end);
            #endif
        #elif defined MANUAL_SCATTER
            long remainder_loop_start = v_end;
            loop_end = loop_end_orig;
            for (int x=0; x<NVAR; x++) {
                for (long i=remainder_loop_start; i<loop_end; i++) {
                    fluxes_a[x][i-remainder_loop_start] = 0.0;
                    fluxes_b[x][i-remainder_loop_start] = 0.0;
                    #ifdef MANUAL_GATHER
                        variables_a[x][i-remainder_loop_start] = variables[edges[i].a*NVAR+x];
                        variables_b[x][i-remainder_loop_start] = variables[edges[i].b*NVAR+x];
                    #endif
                }
            }
            #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
                // Each thread already has its own remainders to process, no need to 
                // distribute workload
            #endif
        #endif

        #pragma omp simd safelen(1)
        for (long i=remainder_loop_start; i<loop_end; i++)
        {
            compute_flux_edge_kernel_crippled(
                #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                    edge_weights[i],
                #else
                    edges[i].x, 
                    edges[i].y, 
                    edges[i].z,
                #endif
                #if defined SIMD && defined MANUAL_GATHER
                    variables_a, 
                    variables_b,
                #else
                    &variables[edges[i].a*NVAR],
                    &variables[edges[i].b*NVAR],
                #endif
                #ifdef FLUX_FISSION
                    &edge_variables[i*NVAR]
                #else
                    #ifdef MANUAL_SCATTER
                        i-remainder_loop_start,
                        fluxes_a, 
                        fluxes_b
                    #else
                        &fluxes[edges[i].a*NVAR],
                        &fluxes[edges[i].b*NVAR]
                    #endif
                #endif
                );
        }
        #ifdef MANUAL_SCATTER
            // Write out fluxes:
            for (long i=remainder_loop_start; i<loop_end; i++) {
                for (int x=0; x<NVAR; x++) {
                    fluxes[edges[i].a*NVAR+x] += fluxes_a[x][i-remainder_loop_start];
                    fluxes[edges[i].b*NVAR+x] += fluxes_b[x][i-remainder_loop_start];
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
