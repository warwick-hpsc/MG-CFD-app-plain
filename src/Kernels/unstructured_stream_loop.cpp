#include "unstructured_stream_loop.h"
#include "unstructured_stream_kernel.h"

#include "cfd_loops.h"

#include "papi_funcs.h"
#include "timer.h"
#include "loop_stats.h"

// Performs same data movement as compute_flux_edge() but with minimal arithmetic. 
// Thus measures memory-bound of compute_flux_edge().
void unstructured_stream_loop(
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
    #endif
    )
{
    log("Performing unstructured_stream_loop()");
    current_kernel = UNSTRUCTURED_STREAM;

    long loop_start = first_edge;
    long loop_end = loop_start + nedges;

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

    #ifdef __clang__
        #pragma clang loop vectorize(disable)
    #else
        #pragma omp simd safelen(1)
    #endif
    #pragma nounroll
    for (long i=loop_start; i<loop_end; i++)
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

    #ifdef TIME
    stop_timer();
    #endif
    #ifdef PAPI
    stop_papi();
    #endif

    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        }
    #endif

    log("unstructured_stream_loop() complete");
}
