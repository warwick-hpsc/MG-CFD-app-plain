#include "indirect_rw_kernel.h"
#include "kernels.h"

#include "papi_funcs.h"
#include "timer.h"
#include "loop_stats.h"

// Indirect R/W kernel
// - performs same data movement as compute_flux_edge() but with minimal arithmetic, 
//   measuring upper bound on performance achievable by compute_flux_edge()
void indirect_rw(
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
    #endif
    )
{
    log("Performing indirect RW");
    current_kernel = INDIRECT_RW;

    int loop_start = first_edge;
    int loop_end = loop_start + nedges;

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
            #endif
        #endif
    #endif
    for (int i=loop_start; i<loop_end; i++)
    {
        #include "indirect_rw_kernel.elemfunc.c"
    }
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

    log("Indirect RW complete");
}
