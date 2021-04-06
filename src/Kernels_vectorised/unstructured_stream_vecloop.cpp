#include "unstructured_stream_vecloop.h"
#include "unstructured_stream_veckernel.h"
#include "unstructured_stream_kernel.h"
#include "cfd_loops.h"

#include "papi_funcs.h"
#include "timer.h"
#include "loop_stats.h"

// Performs same data movement as compute_flux_edge() but with minimal arithmetic. 
// Thus measures memory-bound of compute_flux_edge().
void unstructured_stream_vecloop(
    long first_edge,
    long nedges,
    const long *restrict edge_nodes, 
    const double *restrict edge_vectors,
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
    log("Performing unstructured_stream_vecloop()");
    current_kernel = UNSTRUCTURED_STREAM;

    long loop_start = first_edge;
    long loop_end = loop_start + nedges;

    #if defined COLOURED_CONFLICT_AVOIDANCE && (!defined FLUX_FISSION)
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
    #endif

    #ifdef PAPI
    start_papi();
    #endif
    #ifdef TIME
    start_timer();
    #endif
    record_iters(loop_start, loop_end);

    #ifdef FLUX_FISSION
        // SIMD is safe
        #ifdef DBLS_PER_SIMD
            SIMD_LOOP(DBLS_PER_SIMD)
        #else
            SIMD_LOOP_AUTO
        #endif
    #else
        // Conflict avoidance is required for safe SIMD
        #if defined COLOURED_CONFLICT_AVOIDANCE
            SIMD_LOOP(DBLS_PER_SIMD)
        #elif defined MANUAL_SCATTER
            const long loop_end_orig = loop_end;
            long v_start = loop_start;
            long v_end = loop_start + ((loop_end-loop_start)/DBLS_PER_SIMD)*DBLS_PER_SIMD;
            for (long v=v_start; v<v_end; v+=DBLS_PER_SIMD) {
                #ifdef MANUAL_GATHER
                    for (int n=0; n<DBLS_PER_SIMD; n++) {
                        UNROLL_LOOP(NVAR)
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
                loop_start = 0;
                loop_end = DBLS_PER_SIMD;

                SIMD_LOOP(DBLS_PER_SIMD)

        #elif defined USE_AVX512CD
            // Always prefer using OMP pragma to vectorise, gives better performance 
            // than default auto-vectoriser triggered by absent pragma
            #pragma omp simd simdlen(DBLS_PER_SIMD)

        #endif
    #endif
    for (long i=loop_start; i<loop_end; i++)
    {
        #if defined MANUAL_GATHER || defined MANUAL_SCATTER
            const long edge_idx = v+i;
            const int simd_idx = i;
        #else
            const long edge_idx = i;
        #endif
        #if defined USE_AVX512CD
            // For Intel AVX-512-CD auto-vectorizer to act, I need to 
            // directly include the kernel source here rather than 
            // call an inlined kernel function.
            #include "unstructured_stream_veckernel.elemfunc.c"
        #else
            unstructured_stream_veckernel(
                #if defined MANUAL_GATHER || defined MANUAL_SCATTER
                    simd_idx,
                #endif

                #if defined MANUAL_GATHER
                    #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                        simd_ewt,
                    #endif
                    simd_edge_vectors,
                    simd_variables_a, 
                    simd_variables_b,
                #else
                    #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                        edge_weights[edge_idx],
                    #endif
                    edge_vectors[edge_idx*NDIM], edge_vectors[edge_idx*NDIM+1], edge_vectors[edge_idx*NDIM+2],
                    &variables[edge_nodes[edge_idx*2]  *NVAR],
                    &variables[edge_nodes[edge_idx*2+1]*NVAR],
                #endif

                #if defined MANUAL_SCATTER
                    simd_fluxes_a, 
                    simd_fluxes_b
                #elif defined FLUX_FISSION
                    &edge_variables[edge_idx*NVAR]
                #else
                    &fluxes[edge_nodes[edge_idx*2  ]*NVAR],
                    &fluxes[edge_nodes[edge_idx*2+1]*NVAR]
                #endif
                );
        #endif
    }

    #if defined MANUAL_SCATTER && (!defined FLUX_FISSION)
        // Write out fluxes:
            for (long n=0; n<DBLS_PER_SIMD; n++) {
                UNROLL_LOOP(NVAR)
                for (long x=0; x<NVAR; x++) {
                    fluxes[edge_nodes[(v+n)*2]  *NVAR+x] += simd_fluxes_a[x][n];
                    fluxes[edge_nodes[(v+n)*2+1]*NVAR+x] += simd_fluxes_b[x][n];
                }
            }
        } // Close outer loop over SIMD blocks
    #endif

    #if (defined COLOURED_CONFLICT_AVOIDANCE || defined MANUAL_SCATTER) && (!defined FLUX_FISSION)
        // Compute fluxes of 'remainder' edges without SIMD:
        #ifdef COLOURED_CONFLICT_AVOIDANCE
            long remainder_loop_start = serial_section_start;
            loop_end = first_edge + nedges;
        #elif defined MANUAL_SCATTER
            long remainder_loop_start = v_end;
            loop_end = loop_end_orig;
        #endif
        NOSIMD
        for (long i=remainder_loop_start; i<loop_end; i++)
        {
            unstructured_stream_kernel(
                #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                    edge_weights[i],
                #endif
                edge_vectors[i*NDIM], edge_vectors[i*NDIM+1], edge_vectors[i*NDIM+2],
                &variables[edge_nodes[i*2]  *NVAR],
                &variables[edge_nodes[i*2+1]*NVAR],
                #ifdef FLUX_FISSION
                    &edge_variables[i*NVAR]
                #else
                    &fluxes[edge_nodes[i*2  ]*NVAR],
                    &fluxes[edge_nodes[i*2+1]*NVAR]
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

    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        }
    #endif

    log("unstructured_stream_vecloop() complete");
}
