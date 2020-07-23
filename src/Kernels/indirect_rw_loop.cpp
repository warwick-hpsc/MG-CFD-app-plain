#include "indirect_rw_loop.h"
#include "indirect_rw_kernel.h"
#include "cfd_loops.h"

#include "papi_funcs.h"
#include "timer.h"
#include "loop_stats.h"

// Indirect R/W kernel
// - performs same data movement as compute_flux_edge() but with minimal arithmetic. 
//   Measures upper bound on performance achievable by compute_flux_edge()
void indirect_rw(
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
    log("Performing indirect RW");
    current_kernel = INDIRECT_RW;

    long loop_start = first_edge;
    long loop_end = loop_start + nedges;

    // edges = (edge_neighbour*)__builtin_assume_aligned(edges, 64);

    #if defined SIMD && defined COLOURED_CONFLICT_AVOIDANCE && (!defined FLUX_FISSION)
        loop_end = serial_section_start;
    #endif

    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        #pragma omp parallel firstprivate(loop_start, loop_end)
        {
            openmp_distribute_loop_iterations(&loop_start, &loop_end);
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

    #ifdef PAPI
    start_papi();
    #endif
    #ifdef TIME
    start_timer();
    #endif
    record_iters(loop_start, loop_end);

    #ifndef SIMD
        #ifdef __clang__
            #pragma clang loop vectorize(disable)
        #else
            #pragma omp simd safelen(1)
        #endif
    #else
        #ifdef FLUX_FISSION
            // SIMD is safe
            #ifdef __clang__
                #pragma clang loop vectorize_width(DBLS_PER_SIMD)
            #else
                #pragma omp simd simdlen(DBLS_PER_SIMD)
            #endif
        #else
            // Conflict avoidance is required for safe SIMD
            #if defined COLOURED_CONFLICT_AVOIDANCE
                #ifdef __clang__
                    #pragma clang loop vectorize_width(DBLS_PER_SIMD)
                    #pragma nounroll
                #else
                    #pragma omp simd simdlen(DBLS_PER_SIMD)
                #endif
            #elif defined MANUAL_SCATTER
                const long loop_end_orig = loop_end;
                long v_start = loop_start;
                long v_end = loop_start + ((loop_end-loop_start)/MANUAL_WIDTH)*MANUAL_WIDTH;
                for (long v=v_start; v<v_end; v+=MANUAL_WIDTH) {
                    #ifdef MANUAL_GATHER
                        #ifdef __clang__
                            // Preventing unrolling of outer loop helps assembly-loop-extractor
                            #pragma nounroll
                            #pragma clang loop vectorize(disable)
                        #else
                            #pragma omp simd safelen(1)
                        #endif
                        for (int n=0; n<MANUAL_WIDTH; n++) {
                            #ifdef __clang__
                                #pragma unroll
                            #endif
                            for (int x=0; x<NVAR; x++) {
                                variables_a[x][n] = variables[edges[v+n].a*NVAR+x];
                                variables_b[x][n] = variables[edges[v+n].b*NVAR+x];
                            }
                        }
                    #endif
                    loop_start = v;
                    loop_end = v+MANUAL_WIDTH;

                    #ifdef __clang__
                        // __builtin_assume(flux_loop_start%DBLS_PER_SIMD == 0);
                        // __builtin_assume(loop_end%DBLS_PER_SIMD == 0);
                        #pragma clang loop vectorize_width(DBLS_PER_SIMD)
                    #else
                        #pragma omp simd simdlen(DBLS_PER_SIMD)
                    #endif

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
            #include "indirect_rw_kernel.elemfunc.c"
        #elif defined __clang__
            // Clang is not inlining call to indirect_rw_kernel(), 
            // so manually inline here:
            #include "indirect_rw_kernel.elemfunc.c"
        #elif defined __GNUC__
            // Inlining call makes GNU vector log easier to parse
            #include "indirect_rw_kernel.elemfunc.c"
        #else
            indirect_rw_kernel(
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

    #if defined SIMD && defined MANUAL_SCATTER && (!defined FLUX_FISSION)
        // Write out fluxes:
            #ifdef __clang__
                // Preventing unrolling of outer loop helps assembly-loop-extractor
                #pragma nounroll
            #endif
            for (long n=0; n<MANUAL_WIDTH; n++) {
                #ifdef __clang__
                    #pragma clang loop vectorize(disable)
                    #pragma unroll
                #else
                    #pragma omp simd safelen(1)
                #endif
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

        #ifdef __clang__
            #pragma clang loop vectorize(disable)
            #pragma nounroll
        #else
            #pragma omp simd safelen(1)
        #endif
        for (long i=remainder_loop_start; i<loop_end; i++)
        {
            indirect_rw_kernel(
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
                for (int v=0; v<NVAR; v++) {
                    fluxes[edges[i].a*NVAR+v] += fluxes_a[v][i-remainder_loop_start];
                    fluxes[edges[i].b*NVAR+v] += fluxes_b[v][i-remainder_loop_start];
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

    log("Indirect RW complete");
}
