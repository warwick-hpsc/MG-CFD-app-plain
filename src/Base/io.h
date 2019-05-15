#ifndef IO_H
#define IO_H

#include "common.h"

void read_grid(
    const char* data_file_name, 
    int* nel, 
    double** volumes, 
    int* number_of_edges, 
    int* num_internal_edges, 
    int* num_boundary_edges, 
    int* num_wall_edges, 
    int* internal_edges_start, 
    int* boundary_edges_start, 
    int* wall_edges_start, 
    edge_neighbour** edges, 
    double3** coords);

void dump(
    const double *restrict variables, 
    int nel, int level);

#endif

