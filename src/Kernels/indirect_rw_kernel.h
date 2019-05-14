#ifndef INDIRECT_RW_H
#define INDIRECT_RW_H

#include "common.h"

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
    );

#endif
