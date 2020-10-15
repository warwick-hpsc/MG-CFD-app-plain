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
    const long *restrict edge_nodes, 
    #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
        const double *restrict edge_weights,
    #endif
    const double *restrict edge_vectors,
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

    #if defined SIMD && defined COLOURED_CONFLICT_AVOIDANCE && (!defined FLUX_FISSION)
        loop_end = serial_section_start;
    #endif

    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        #pragma omp parallel firstprivate(loop_start, loop_end)
        {
            openmp_distribute_loop_iterations(&loop_start, &loop_end);
    #endif

    #if defined MANUAL_SCATTER
        double simd_fluxes_a[NVAR][DBLS_PER_SIMD];
        double simd_fluxes_b[NVAR][DBLS_PER_SIMD];
        #ifdef MANUAL_GATHER
            double simd_variables_a[NVAR][DBLS_PER_SIMD];
            double simd_variables_b[NVAR][DBLS_PER_SIMD];
        #endif
        double simd_edge_vectors[NDIM][DBLS_PER_SIMD];
        double simd_ewt[DBLS_PER_SIMD];

        for (int n=0; n<DBLS_PER_SIMD; n++) {
            for (int x=0; x<NVAR; x++) {
                simd_fluxes_a[x][n] = 0.0;
                simd_fluxes_b[x][n] = 0.0;
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
        #pragma nounroll
    #else
        #ifdef FLUX_FISSION
            // SIMD is safe
            #ifdef __clang__
                #pragma clang loop vectorize_width(DBLS_PER_SIMD)
            #else
                #pragma omp simd simdlen(DBLS_PER_SIMD)
            #endif
            #pragma nounroll
        #else
            // Conflict avoidance is required for safe SIMD
            #if defined COLOURED_CONFLICT_AVOIDANCE
                #ifdef __clang__
                    #pragma clang loop vectorize_width(DBLS_PER_SIMD) interleave(disable)
                #else
                    #pragma omp simd simdlen(DBLS_PER_SIMD)
                #endif
            #elif defined MANUAL_SCATTER
                const long loop_end_orig = loop_end;
                long v_start = loop_start;
                long v_end = loop_start + ((loop_end-loop_start)/DBLS_PER_SIMD)*DBLS_PER_SIMD;
                for (long v=v_start; v<v_end; v+=DBLS_PER_SIMD) {
                    #ifdef MANUAL_GATHER
                        #ifdef __clang__
                            // Warning: if you ask Clang to vectorize this loop containing 
                            //          indirect reads, it will do it BUT the resulting 
                            //          assembly sequence is 10x longer than non-SIMD - 
                            //          a mangled mess of shuffles and moves, serial and packed.
                            //          Doesn't affect performance, but DOES affect ability of 
                            //          assembly-loop-extractor to identify this loop.
                            //
                            //          Possibly only Intel can handle this loop nicely. 
                            //
                            // #pragma clang loop vectorize_width(DBLS_PER_SIMD) interleave(disable)
                        #else
                            #pragma omp simd simdlen(DBLS_PER_SIMD)
                        #endif
                        for (int n=0; n<DBLS_PER_SIMD; n++) {
                            #pragma unroll
                            for (int x=0; x<NVAR; x++) {
                                simd_variables_a[x][n] = variables[edge_nodes[(v+n)*2]  *NVAR+x];
                                simd_variables_b[x][n] = variables[edge_nodes[(v+n)*2+1]*NVAR+x];
                            }

                            #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                                simd_ewt[n] = edge_weights[v+n];
                            #endif
                            simd_edge_vectors[0][n] = edge_vectors[(v+n)*NDIM];
                            simd_edge_vectors[1][n] = edge_vectors[(v+n)*NDIM+1];
                            simd_edge_vectors[2][n] = edge_vectors[(v+n)*NDIM+2];
                        }
                    #endif
                    loop_start = v;
                    loop_end = v+DBLS_PER_SIMD;

                    #pragma omp simd simdlen(DBLS_PER_SIMD)

            #elif defined USE_AVX512CD
                // Always prefer using OMP pragma to vectorise, gives better performance 
                // than default auto-vectoriser triggered by absent pragma
                #pragma omp simd simdlen(DBLS_PER_SIMD)

            #else
                #pragma omp simd safelen(1)
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
        #elif defined __GNUC__ && ! defined __clang__
            // Inlining call makes GNU vector log easier to parse
            #include "indirect_rw_kernel.elemfunc.c"
        #else
            indirect_rw_kernel(
                #if defined SIMD && (defined MANUAL_GATHER || defined MANUAL_SCATTER)
                    i-loop_start,
                #endif

                #if defined SIMD && defined MANUAL_GATHER
                    #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                        simd_ewt,
                    #endif
                    simd_edge_vectors,
                    simd_variables_a, 
                    simd_variables_b,
                #else
                    #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                        edge_weights[i],
                    #endif
                    edge_vectors[i*NDIM], edge_vectors[i*NDIM+1], edge_vectors[i*NDIM+2],
                    &variables[edge_nodes[i*2]  *NVAR],
                    &variables[edge_nodes[i*2+1]*NVAR],
                #endif

                #if defined SIMD && defined MANUAL_SCATTER
                    simd_fluxes_a, 
                    simd_fluxes_b
                #elif defined FLUX_FISSION
                    &edge_variables[i*NVAR]
                #else
                    &fluxes[edge_nodes[i*2  ]*NVAR],
                    &fluxes[edge_nodes[i*2+1]*NVAR]
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
                    fluxes[edge_nodes[(v+n)*2]  *NVAR+x] += simd_fluxes_a[x][n];
                    fluxes[edge_nodes[(v+n)*2+1]*NVAR+x] += simd_fluxes_b[x][n];
                    simd_fluxes_a[x][n] = 0.0;
                    simd_fluxes_b[x][n] = 0.0;
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
                    simd_fluxes_a[x][simd_idx] = 0.0;
                    simd_fluxes_b[x][simd_idx] = 0.0;
                    #ifdef MANUAL_GATHER
                        simd_variables_a[x][simd_idx] = variables[edge_nodes[i*2  ]*NVAR+x];
                        simd_variables_b[x][simd_idx] = variables[edge_nodes[i*2+1]*NVAR+x];
                    #endif
                }
                #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                    simd_ewt[simd_idx] = edge_weights[i];
                #endif
                simd_edge_vectors[0][simd_idx] = edge_vectors[i*NDIM];
                simd_edge_vectors[1][simd_idx] = edge_vectors[i*NDIM+1];
                simd_edge_vectors[2][simd_idx] = edge_vectors[i*NDIM+2];
            }
            #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
                // Each thread already has its own remainders to process, no need to 
                // distribute workload
            #endif
        #endif

        #pragma omp simd safelen(1)
        for (long i=remainder_loop_start; i<loop_end; i++)
        {
            indirect_rw_kernel(
                #if defined SIMD && (defined MANUAL_GATHER || defined MANUAL_SCATTER)
                    i-remainder_loop_start,
                #endif

                #if defined SIMD && defined MANUAL_GATHER
                    #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                        simd_ewt,
                    #endif
                    simd_edge_vectors,
                    simd_variables_a, 
                    simd_variables_b,
                #else
                    #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                        edge_weights[i],
                    #endif
                    edge_vectors[i*NDIM], edge_vectors[i*NDIM+1], edge_vectors[i*NDIM+2],
                    &variables[edge_nodes[i*2  ]*NVAR],
                    &variables[edge_nodes[i*2+1]*NVAR],
                #endif

                #if defined SIMD && defined MANUAL_SCATTER
                    simd_fluxes_a, 
                    simd_fluxes_b
                #elif defined FLUX_FISSION
                    &edge_variables[i*NVAR]
                #else
                    &fluxes[edge_nodes[i*2  ]*NVAR],
                    &fluxes[edge_nodes[i*2+1]*NVAR]
                #endif
                );
        }
        #ifdef MANUAL_SCATTER
            // Write out fluxes:
            for (long i=remainder_loop_start; i<loop_end; i++) {
                for (int v=0; v<NVAR; v++) {
                    fluxes[edge_nodes[i*2  ]*NVAR+v] += simd_fluxes_a[v][i-remainder_loop_start];
                    fluxes[edge_nodes[i*2+1]*NVAR+v] += simd_fluxes_b[v][i-remainder_loop_start];
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
