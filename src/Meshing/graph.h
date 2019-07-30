#ifndef GRAPH_H
#define GRAPH_H

#include "common.h"

// Count degree of each point.
int* count_point_degrees(
	int npoints, 
	const edge_neighbour* edges, 
	int edges_start_idx, 
	int edges_end_idx);

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

// Identify edges connected to each point.
int** construct_point_to_edge_mapping(
	const int* point_degree, 
	int npoints, 
	const edge_neighbour* edges, 
	int edges_start_idx, 
	int edges_end_idx);

// Identify edgewise neighbours
int* count_edge_degrees(
	const int* point_degree,
	int** point_to_edges, 
	const edge_neighbour* edges, 
	int nedges, 
	int edges_start, 
	int edges_end);

// Generate mapping of edge_id -> {edge_ids}
int** identify_edgewise_neighbours(
	int** point_to_edges, 
	const int* point_degree,
	const int* edge_degree,
	const edge_neighbour* edges,
	int nedges,
	int edges_start_idx,
	int edges_end_idx);

#endif