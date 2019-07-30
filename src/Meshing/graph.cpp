#include <stdio.h>
#include <stdlib.h>

#include "graph.h"

// Count internal-degree of each point.
int* count_point_internal_degrees(
	int npoints, 
	const edge_neighbour* edges, 
	int edges_start_idx, 
	int edges_end_idx)
{
	int* point_degree = alloc<int>(npoints);
	for (int i=0; i<npoints; i++) point_degree[i] = 0;
	for (int e=edges_start_idx; e<=edges_end_idx; e++) {
		if (edges[e].a >= 0 && edges[e].b >= 0) {
			point_degree[edges[e].a]++;
			point_degree[edges[e].b]++;
		}
	}

	return point_degree;
}

// Identify point-wise internal neighbours
int** identify_pointwise_internal_neighbours(
	const int* point_degree, 
	int npoints, 
	const edge_neighbour* edges, 
	int edges_start_idx, 
	int edges_end_idx)
{
	int** point_neighbours = alloc<int*>(npoints);
	for (int i=0; i<npoints; i++) {
		if (point_degree[i] == 0) {
			point_neighbours[i] = NULL;
		} else {
			point_neighbours[i] = alloc<int>(point_degree[i]);
		}
	}

	int a, b;
	int* ptrs = alloc<int>(npoints);
	for (int i=0; i<npoints; i++) ptrs[i] = 0;
	for (int e=edges_start_idx; e<=edges_end_idx; e++) {
		a = edges[e].a;
		b = edges[e].b;
		if (a >= 0 && b >= 0) {
			point_neighbours[a][ptrs[a]] = b;
			ptrs[a]++;
			point_neighbours[b][ptrs[b]] = a;
			ptrs[b]++;
		}
	}

	dealloc<int>(ptrs);
	return point_neighbours;
}
