#ifndef FLUX_KERNELS_H
#define FLUX_KERNELS_H

#include "common.h"

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
    );

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
    );

void compute_flux_edge_loop(
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
    );

#ifdef FLUX_CRIPPLE
void compute_flux_edge_crippled_loop(
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
    );
#endif

#endif
