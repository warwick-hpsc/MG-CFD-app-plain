#ifndef INDIRECT_RW_VECLOOP_H
#define INDIRECT_RW_VECLOOP_H

#include "common.h"

void indirect_rw_vecloop(
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
    );

#endif
