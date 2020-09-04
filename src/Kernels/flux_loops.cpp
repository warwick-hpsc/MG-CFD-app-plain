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

    #if defined SIMD && defined COLOURED_CONFLICT_AVOIDANCE && (!defined FLUX_FISSION)
        loop_end = serial_section_start;
    #endif

    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        #pragma omp parallel firstprivate(flux_loop_start, loop_end)
        {
            openmp_distribute_loop_iterations(&flux_loop_start, &loop_end);
    #endif

    #if defined MANUAL_SCATTER
        double fluxes_a[NVAR][DBLS_PER_SIMD];
        double fluxes_b[NVAR][DBLS_PER_SIMD];
        #ifdef MANUAL_GATHER
            double variables_a[NVAR][DBLS_PER_SIMD];
            double variables_b[NVAR][DBLS_PER_SIMD];
        #endif
        double edge_vectors[NDIM][DBLS_PER_SIMD];

        for (int n=0; n<DBLS_PER_SIMD; n++) {
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
        #ifdef __clang__
            #pragma clang loop vectorize(disable)
        #else
            #pragma omp simd safelen(1)
        #endif
    #else
        #ifdef FLUX_FISSION
            // SIMD is safe
            #pragma omp simd simdlen(DBLS_PER_SIMD)
        #else
            // Conflict avoidance is required for safe SIMD
            #if defined COLOURED_CONFLICT_AVOIDANCE
                #pragma omp simd simdlen(DBLS_PER_SIMD)
            #elif defined MANUAL_SCATTER
                const long loop_end_orig = loop_end;
                long v_start = flux_loop_start;
                long v_end = flux_loop_start + ((loop_end-flux_loop_start)/DBLS_PER_SIMD)*DBLS_PER_SIMD;
                for (long v=v_start; v<v_end; v+=DBLS_PER_SIMD) {
                    #ifdef MANUAL_GATHER
                        #ifdef __clang__
                            // When Clang encouters "#pragma omp simd" it turns on its loop unroller, 
                            // even when '-fno-unroll-loops' flag was passed to turn it off.
                            // Frustrating as this gather loop needs to be vectorised. 
                            // So use Clang-specific pragma instead:
                            #pragma clang loop vectorize_width(DBLS_PER_SIMD)
                        #else
                            #pragma omp simd simdlen(DBLS_PER_SIMD)
                        #endif
                        // Preventing unrolling of outer loop helps assembly-loop-extractor
                        #pragma nounroll
                        for (int n=0; n<DBLS_PER_SIMD; n++) {
                            #pragma unroll
                            for (int x=0; x<NVAR; x++) {
                                variables_a[x][n] = variables[edges[v+n].a*NVAR+x];
                                variables_b[x][n] = variables[edges[v+n].b*NVAR+x];
                            }
                            edge_vectors[0][n] = edges[v+n].x;
                            edge_vectors[1][n] = edges[v+n].y;
                            edge_vectors[2][n] = edges[v+n].z;
                        }
                    #endif
                    flux_loop_start = v;
                    loop_end = v+DBLS_PER_SIMD;

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
        #elif defined __GNUC__ && ! defined __clang__
            // Inlining call makes GNU vector log easier to parse
            #include "flux_kernel.elemfunc.c"
        #else
            compute_flux_edge_kernel(
                #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                    edge_weights[i],
                #endif
                #if defined SIMD && defined MANUAL_GATHER
                    edge_vectors,
                    variables_a, 
                    variables_b,
                #else
                    edges[i].x, edges[i].y, edges[i].z,
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
            for (long n=0; n<DBLS_PER_SIMD; n++) {
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
            for (long i=remainder_loop_start; i<loop_end; i++) {
                const int simd_idx = i-remainder_loop_start;
                for (int x=0; x<NVAR; x++) {
                    fluxes_a[x][simd_idx] = 0.0;
                    fluxes_b[x][simd_idx] = 0.0;
                    #ifdef MANUAL_GATHER
                        variables_a[x][simd_idx] = variables[edges[i].a*NVAR+x];
                        variables_b[x][simd_idx] = variables[edges[i].b*NVAR+x];
                    #endif
                }
                edge_vectors[0][simd_idx] = edges[i].x;
                edge_vectors[1][simd_idx] = edges[i].y;
                edge_vectors[2][simd_idx] = edges[i].z;
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
                #if defined SIMD && defined MANUAL_GATHER
                    edge_vectors,
                    variables_a, 
                    variables_b,
                #else
                    edges[i].x, edges[i].y, edges[i].z,
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
        double fluxes_a[NVAR][DBLS_PER_SIMD];
        double fluxes_b[NVAR][DBLS_PER_SIMD];
        #ifdef MANUAL_GATHER
            double variables_a[NVAR][DBLS_PER_SIMD];
            double variables_b[NVAR][DBLS_PER_SIMD];
        #endif
        for (int n=0; n<DBLS_PER_SIMD; n++) {
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
                long v_end = flux_loop_start + ((loop_end-flux_loop_start)/DBLS_PER_SIMD)*DBLS_PER_SIMD;
                #pragma nounroll
                for (long v=v_start; v<v_end; v+=DBLS_PER_SIMD) {
                    #ifdef MANUAL_GATHER
                        // Preventing unrolling of outer loop helps assembly-loop-extractor
                        #pragma nounroll
                        #pragma omp simd safelen(1)
                        for (int n=0; n<DBLS_PER_SIMD; n++) {
                            #pragma unroll
                            for (int x=0; x<NVAR; x++) {
                                variables_a[x][n] = variables[edges[v+n].a*NVAR+x];
                                variables_b[x][n] = variables[edges[v+n].b*NVAR+x];
                            }
                        }
                    #endif
                    flux_loop_start = v;
                    loop_end = v+DBLS_PER_SIMD;

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
            for (long n=0; n<DBLS_PER_SIMD; n++) {
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
