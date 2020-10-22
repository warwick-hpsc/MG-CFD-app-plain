#ifndef IO_H
#define IO_H

#include "common.h"

void read_grid(
    const char* data_file_name, 
    long* nel, 
    double** volumes, 
    long* number_of_edges, 
    long* num_internal_edges, 
    long* num_boundary_edges, 
    long* num_wall_edges, 
    long* internal_edges_start, 
    long* boundary_edges_start, 
    long* wall_edges_start, 
    edge_neighbour** edges, 
    double3** coords);

void dump(
    const double *restrict variables, 
    long nel, int level);

#endif

