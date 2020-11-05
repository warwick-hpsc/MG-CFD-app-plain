#include "kernel_wrappers.h"

#include "flux_loops.h"
#include "flux_vecloops.h"

#include "unstructured_compute_loop.h"
#include "unstructured_compute_vecloop.h"

#include "unstructured_stream_loop.h"
#include "unstructured_stream_vecloop.h"

#include "compute_stream_loop.h"

void compute_flux_edge(
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
    #ifdef SIMD
        compute_flux_edge_vecloop(
            first_edge,
            nedges,
            edge_nodes, 
            edge_vectors,
            #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                edge_weights,
            #endif
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
            edge_vectors,
            #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                edge_weights,
            #endif
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

void unstructured_compute(
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
    #ifdef SIMD
        unstructured_compute_vecloop(
            first_edge,
            nedges,
            edge_nodes, 
            edge_vectors,
            #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                edge_weights,
            #endif
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
        unstructured_compute_loop(
            first_edge,
            nedges,
            edge_nodes, 
            edge_vectors,
            #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                edge_weights,
            #endif
            variables, 
            #ifdef FLUX_FISSION
                edge_variables
            #else
                fluxes
            #endif
            );
    #endif
}

void unstructured_stream(
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
    #ifdef SIMD
        unstructured_stream_vecloop(
            first_edge,
            nedges,
            edge_nodes, 
            edge_vectors,
            #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                edge_weights,
            #endif
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
        unstructured_stream_loop(
            first_edge,
            nedges,
            edge_nodes, 
            edge_vectors,
            #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                edge_weights,
            #endif
            variables, 
            #ifdef FLUX_FISSION
                edge_variables
            #else
                fluxes
            #endif
            );
    #endif
}

void compute_stream(
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
    // #ifdef SIMD
    //     compute_stream_vecloop(
    //         first_edge,
    //         nedges,
    //         edge_nodes, 
    //         edge_vectors,
    //         #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
    //             edge_weights,
    //         #endif
    //         variables, 
    //         #ifdef FLUX_FISSION
    //             edge_variables
    //         #else
    //             fluxes
    //             #ifdef COLOURED_CONFLICT_AVOIDANCE
    //             , serial_section_start
    //             #endif
    //         #endif
    //         );
    // #else
        compute_stream_loop(
            first_edge,
            nedges,
            edge_nodes, 
            edge_vectors,
            #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                edge_weights,
            #endif
            variables, 
            #ifdef FLUX_FISSION
                edge_variables
            #else
                fluxes
            #endif
            );
    // #endif
}
