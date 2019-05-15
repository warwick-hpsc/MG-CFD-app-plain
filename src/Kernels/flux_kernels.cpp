#include "flux_kernels.h"
#include "kernels.h"

#include "papi_funcs.h"
#include "timer.h"
#include "loop_stats.h"

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
        #include "flux_boundary.elemfunc.c"
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
        #include "flux_wall.elemfunc.c"
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
    #endif
    )
{
    log("Computing internal flux");
    current_kernel = COMPUTE_FLUX_EDGE;

    int loop_start = first_edge;
    int loop_end = loop_start + nedges;

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
                // TODO: Insert the following pragma into flux.elemfunc.c to 
                //       enable safe AVX-512-CD SIMD:
                // #pragma omp ordered simd overlap(...)
            #endif
        #endif
    #endif
    for (int i=loop_start; i<loop_end; i++)
    {
        #include "flux.elemfunc.c"
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

    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        }
    #endif

    log("Internal flux compute complete");
}

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
    #endif
    )
{
    log("Computing internal flux crippled");
    current_kernel = COMPUTE_FLUX_EDGE;

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
        #include "flux-crippled.elemfunc.c"
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
}
