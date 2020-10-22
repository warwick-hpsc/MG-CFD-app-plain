#include "flux_loops.h"
#include "cfd_loops.h"

#include "flux_kernel.h"
#ifdef FLUX_CRIPPLE
    #include "flux_kernel_crippled.h"
#endif
#include "flux_boundary_kernel.h"
#include "flux_wall_kernel.h"

#include "papi_funcs.h"
#include "timer.h"
#include "loop_stats.h"

void compute_flux_edge_loop(
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
    log("Computing internal flux");
    current_kernel = COMPUTE_FLUX_EDGE;

    long flux_loop_start = first_edge;
    long loop_end = flux_loop_start + nedges;

    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        #pragma omp parallel firstprivate(flux_loop_start, loop_end)
        {
            openmp_distribute_loop_iterations(&flux_loop_start, &loop_end);
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

    #pragma omp simd safelen(1)
    for (long i=flux_loop_start; i<loop_end; i++)
    {
        compute_flux_edge_kernel(
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
void compute_flux_edge_crippled_loop(
    long first_edge,
    long nedges,
    const long *restrict edge_nodes, 
    #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
        const double *restrict edge_weights,
    #else
        const double *restrict edge_vectors,
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

    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        #pragma omp parallel firstprivate(flux_loop_start, loop_end)
        {
            openmp_distribute_loop_iterations(&flux_loop_start, &loop_end);
    #endif

    #ifdef PAPI
    start_papi();
    #endif
    #ifdef TIME
    start_timer();
    #endif
    record_iters(flux_loop_start, loop_end);

    #pragma omp simd safelen(1)
    for (long i=flux_loop_start; i<loop_end; i++)
    {
        compute_flux_edge_kernel_crippled(
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
}
#endif

void compute_boundary_flux_edge_loop(
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
    #pragma omp simd safelen(1)
    for (long i=loop_start; i<loop_end; i++)
    {
        compute_boundary_flux_edge_kernel(
            edges[i].x, edges[i].y, edges[i].z,
            &variables[edges[i].b*NVAR],
            #ifdef FLUX_FISSION
                &edge_variables[i*NVAR]
            #else
                &fluxes[edges[i].b*NVAR]
            #endif
            );
    }
    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        }
    #endif

    log("Boundary flux compute complete");
}

void compute_wall_flux_edge_loop(
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
    #pragma omp simd safelen(1)
    for (long i=loop_start; i<loop_end; i++)
    {
        compute_wall_flux_edge_kernel(
            edges[i].x, edges[i].y, edges[i].z,
            &variables[edges[i].b*NVAR],
            #ifdef FLUX_FISSION
                &edge_variables[i*NVAR]
            #else
                &fluxes[edges[i].b*NVAR]
            #endif
            );
    }
    #if defined OMP && (defined FLUX_FISSION || defined OMP_SCATTERS)
        }
    #endif

    log("Wall flux compute complete");
}
