#include "kernel_wrappers.h"

#include "flux_vecloops.h"
#include "flux_loops.h"

#include "indirect_rw_loop.h"
#include "indirect_rw_vecloop.h"

void compute_flux_edge(
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
    #ifdef SIMD
        compute_flux_edge_vecloop(
            first_edge,
            nedges,
            edge_nodes, 
            #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                edge_weights,
            #endif
            edge_vectors,
            variables, 
            #ifdef FLUX_FISSION
                edge_variables
            #else
                fluxes
                #ifdef COLOURED_CONFLICT_AVOIDANCE
                , serial_section_start
                #endif
            #endif
            );
    #else
        compute_flux_edge_loop(
            first_edge,
            nedges,
            edge_nodes, 
            #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                edge_weights,
            #endif
            edge_vectors,
            variables, 
            #ifdef FLUX_FISSION
                edge_variables
            #else
                fluxes
            #endif
            );
    #endif
}

void compute_boundary_flux_edge(
    long first_edge,
    long nedges,
    const edge_neighbour *restrict edges, 
    const double *restrict variables_b, 
    #ifdef FLUX_FISSION
        edge *restrict edge_variables
    #else
        double *restrict fluxes_b
    #endif
    )
{
    compute_boundary_flux_edge_loop(
        first_edge,
        nedges,
        edges, 
        variables_b, 
        #ifdef FLUX_FISSION
            edge_variables
        #else
            fluxes_b
        #endif
        );
}

void compute_wall_flux_edge(
    long first_edge,
    long nedges,
    const edge_neighbour *restrict edges, 
    const double *restrict variables_b, 
    #ifdef FLUX_FISSION
        edge *restrict edge_variables
    #else
        double *restrict fluxes_b
    #endif
    )
{
    compute_wall_flux_edge_loop(
        first_edge,
        nedges,
        edges, 
        variables_b, 
        #ifdef FLUX_FISSION
            edge_variables
        #else
            fluxes_b
        #endif
        );
}

#ifdef FLUX_CRIPPLE
void compute_flux_edge_crippled(
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
    #ifdef SIMD
        compute_flux_edge_crippled_vecloop(
            first_edge,
            nedges,
            edge_nodes, 
            #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                edge_weights,
            #endif
            edge_vectors,
            variables, 
            #ifdef FLUX_FISSION
                edge_variables
            #else
                fluxes
                #ifdef COLOURED_CONFLICT_AVOIDANCE
                , serial_section_start
                #endif
            #endif
            );
    #else
        compute_flux_edge_crippled_loop(
            first_edge,
            nedges,
            edge_nodes, 
            #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                edge_weights,
            #endif
            edge_vectors,
            variables, 
            #ifdef FLUX_FISSION
                edge_variables
            #else
                fluxes
            #endif
            );
    #endif
}
#endif


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
    #ifdef SIMD
        indirect_rw_vecloop(
            first_edge,
            nedges,
            edge_nodes, 
            #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                edge_weights,
            #endif
            edge_vectors,
            variables, 
            #ifdef FLUX_FISSION
                edge_variables
            #else
                fluxes
                #ifdef COLOURED_CONFLICT_AVOIDANCE
                , serial_section_start
                #endif
            #endif
            );
    #else
        indirect_rw_loop(
            first_edge,
            nedges,
            edge_nodes, 
            #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                edge_weights,
            #endif
            edge_vectors,
            variables, 
            #ifdef FLUX_FISSION
                edge_variables
            #else
                fluxes
            #endif
            );
    #endif
}
