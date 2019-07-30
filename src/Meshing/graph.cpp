#include <stdio.h>
#include <stdlib.h>

#include "graph.h"

// Count degree of each point.
int* count_point_degrees(
	int npoints, 
	const edge_neighbour* edges, 
	int edges_start_idx, 
	int edges_end_idx)
{
	int* point_degree = alloc<int>(npoints);
	for (int i=0; i<npoints; i++) point_degree[i] = 0;
	for (int e=edges_start_idx; e<=edges_end_idx; e++) {
		if (edges[e].a >= 0) {
			point_degree[edges[e].a]++;
		}
		if (edges[e].b >= 0) {
			point_degree[edges[e].b]++;
		}
	}

	return point_degree;
}

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

// Identify edges connected to each point.
int** construct_point_to_edge_mapping(
	const int* point_degree, 
	int npoints, 
	const edge_neighbour* edges, 
	int edges_start_idx, 
	int edges_end_idx)
{
	int** point_to_edges = alloc<int*>(npoints);
	for (int i=0; i<npoints; i++) {
		if (point_degree[i] == 0)
			point_to_edges[i] = NULL;
		else
			point_to_edges[i] = alloc<int>(point_degree[i]);
	}

	int a, b;
	int* ptrs = alloc<int>(npoints);
	for (int i=0; i<npoints; i++) ptrs[i] = 0;
	for (int e=edges_start_idx; e<=edges_end_idx; e++) {
		a = edges[e].a;
		b = edges[e].b;
		if (a >= 0) {
			point_to_edges[a][ptrs[a]] = e;
			ptrs[a]++;
		}
		if (b >= 0) {
			point_to_edges[b][ptrs[b]] = e;
			ptrs[b]++;
		}
	}

	#ifdef TESTS
		// Check all edges were mapped:
		for (int i=0; i<npoints; i++) {
			if (ptrs[i] != point_degree[i]) {
				printf("ERROR: Only %d of %d edges mapped for point %d\n", ptrs[i], point_degree[i], i);
				printf("%s:%d\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
		}
	#endif

	dealloc<int>(ptrs);

	return point_to_edges;
}

// Identify edgewise neighbours
int* count_edge_degrees(
	const int* point_degree,
	int** point_to_edges, 
	const edge_neighbour* edges, 
	int nedges, 
	int edges_start, 
	int edges_end)
{
	int* edge_degree = alloc<int>(nedges);
	for (int i=0; i<nedges; i++) edge_degree[i] = 0;
	int a, b, en;
	for (int e=edges_start; e<=edges_end; e++) {
		a = edges[e].a;
		if (a >= 0) {
			for (int p=0; p<point_degree[a]; p++) {
				en = point_to_edges[a][p];
				if (en >= edges_start && en <= edges_end) {
					edge_degree[e]++;
				}
			}
			edge_degree[e]--;
			// Subtract 1 to account for edge 'e'.
		}
		b = edges[e].b;
		if (b >= 0) {
			for (int p=0; p<point_degree[b]; p++) {
				en = point_to_edges[b][p];
				if (en >= edges_start && en <= edges_end) {
					edge_degree[e]++;
				}
			}
			edge_degree[e]--;
			// Subtract 1 to account for edge 'e'.
		}
	    #ifdef TESTS
			if (edge_degree[e - edges_start] < 0) {
				printf("ERROR: Have negative edge degree for edge %d\n", e);
				printf("%s:%d\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
	    #endif
	}

	return edge_degree;
}

// Generate mapping of edge_id -> {edge_ids}
int** identify_edgewise_neighbours(
	int** point_to_edges, 
	const int* point_degree,
	const int* edge_degree,
	const edge_neighbour* edges,
	int nedges,
	int edges_start_idx,
	int edges_end_idx)
{
	int** edge_neighbours = alloc<int*>(nedges);
	for (int e=0; e<nedges; e++) {
		if (edge_degree[e] == 0)
			edge_neighbours[e] = NULL;
		else
			edge_neighbours[e] = alloc<int>(edge_degree[e]);
	}

	int* ptrs = alloc<int>(nedges);
	for (int i=0; i<nedges; i++) ptrs[i] = 0;
	int a, b, a2, b2, en, insert_idx;
	for (int e=edges_start_idx; e<=edges_end_idx; e++) {
		a = edges[e].a;
		b = edges[e].b;
		// Identify neighbouring edges at point 'a':
		if (a >= 0) {
			for (int j=0; j<point_degree[a]; j++) {
				en = point_to_edges[a][j];
				if (en == e)
					continue;
				else if (en >= edges_start_idx && 
						 en <= edges_end_idx)
				{
					edge_neighbours[e][ptrs[e]] = en;
					ptrs[e]++;
				    #ifdef TESTS
						if (ptrs[e] > edge_degree[e]) {
							printf("Error: Have identified too many edges as being neighbours of edge %d\n", e);
							printf("%s:%d\n", __FILE__, __LINE__);
							exit(EXIT_FAILURE);
						}
				    #endif
				}
			}
		}
		// Identify neighbouring edges at point 'b':
		if (b >= 0) {
			for (int j=0; j<point_degree[b]; j++) {
				en = point_to_edges[b][j];
				if (en == e)
					continue;
				else if (en >= edges_start_idx && 
						 en <= edges_end_idx)
				{
					edge_neighbours[e][ptrs[e]] = en;
					ptrs[e]++;
				    #ifdef TESTS
						if (ptrs[e] > edge_degree[e]) {
							printf("Error: Have identified too many edges as being neighbours of edge %d\n", e);
							printf("%s:%d\n", __FILE__, __LINE__);
							exit(EXIT_FAILURE);
						}
				    #endif
				}
			}
		}
	}

	#ifdef TESTS
		// Test that all edge-neighbours have been found:
		for (int e=edges_start_idx; e<=edges_end_idx; e++) {
			if (ptrs[e] < edge_degree[e]) {
				printf("Error: Only %d of %d neighbours found for edge %d\n", ptrs[e], edge_degree[e], e);
				printf("%s:%d\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
		}
	#endif

	dealloc<int>(ptrs);
	return edge_neighbours;
}
