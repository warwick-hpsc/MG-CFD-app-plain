#ifndef GRAPH_H
#define GRAPH_H

#include "common.h"

// Count internal-degree of each point.
int* count_point_internal_degrees(
	int npoints, 
	const edge_neighbour* edges, 
	int edges_start_idx, 
	int edges_end_idx);

// Identify point-wise internal neighbours
int** identify_pointwise_internal_neighbours(
	const int* point_degree, 
	int npoints, 
	const edge_neighbour* edges, 
	int edges_start_idx, 
	int edges_end_idx);

#endif