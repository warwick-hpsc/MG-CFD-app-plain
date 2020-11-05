#include "compute_stream_loop.h"
#include "cfd_loops.h"

#include "flux_kernel.h"

#include "papi_funcs.h"
#include "timer.h"
#include "loop_stats.h"

// Performs same cfd compute_flux_edge() but with minimal data movement. 
// Thus measures compute-bound of compute_flux_edge().
void compute_stream_loop(
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
    log("Performing compute_stream_loop()");
    current_kernel = COMPUTE_STREAM;

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

    const int batch = 8;
    #ifdef __clang__
        #pragma clang loop vectorize(disable)
    #else
        #pragma omp simd safelen(1)
    #endif
    #pragma nounroll
    // for (long i=loop_start; i<loop_end; i++)
    for (long i=loop_start; i<loop_end; i+=batch)
    {
        // const long idx = i;
        // Loop over same handful of edges, which should remain in L1 cache:
        // Each edge reads/writes 23 doubles, x8 edges = 1.44 KB
        // const long idx = i%batch;
        for (long idx=loop_start; idx<loop_start+batch; idx++) {
            compute_flux_edge_kernel(
                #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                    edge_weights[idx],
                #endif
                edge_vectors[idx*NDIM], edge_vectors[idx*NDIM+1], edge_vectors[idx*NDIM+2],
                &variables[edge_nodes[idx*2]  *NVAR],
                &variables[edge_nodes[idx*2+1]*NVAR],
                #ifdef FLUX_FISSION
                    &edge_variables[idx*NVAR]
                #else
                    &fluxes[edge_nodes[idx*2  ]*NVAR],
                    &fluxes[edge_nodes[idx*2+1]*NVAR]
                #endif
                );
        }
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

    log("compute_stream_loop() complete");
}
