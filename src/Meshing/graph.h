#ifndef GRAPH_H
#define GRAPH_H

#include "common.h"

// Count degree of each point.
int* count_point_degrees(
	long npoints, 
	const edge_neighbour* edges, 
	long edges_start_idx, 
	long edges_end_idx);

// Count internal-degree of each point.
int* count_point_internal_degrees(
	long npoints, 
	const edge_neighbour* edges, 
	long edges_start_idx, 
	long edges_end_idx);

// Identify point-wise internal neighbours
long** identify_pointwise_internal_neighbours(
	const int* point_degree, 
	long npoints, 
	const edge_neighbour* edges, 
	long edges_start_idx, 
	long edges_end_idx);

// Identify edges connected to each point.
long** construct_point_to_edge_mapping(
	const int* point_degree, 
	long npoints, 
	const edge_neighbour* edges, 
	long edges_start_idx, 
	long edges_end_idx);

// Identify edgewise neighbours
int* count_edge_degrees(
	const int* point_degree,
	long** point_to_edges, 
	const edge_neighbour* edges, 
	long nedges, 
	long edges_start, 
	long edges_end);

// Generate mapping of edge_id -> {edge_ids}
long** identify_edgewise_neighbours(
	long** point_to_edges, 
	const int* point_degree,
	const int* edge_degree,
	const edge_neighbour* edges,
	long nedges,
	long edges_start_idx,
	long edges_end_idx);

#endif

