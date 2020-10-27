#include <stdio.h>
#include <stdlib.h>

#include "reduce_bw.h"
#include "graph.h"
#include "progress.h"

struct tagged_int { int value, tag; };

int compare_two_tagged_ints(const void* a, const void* b) {
	return ((tagged_int*)a)->value - ((tagged_int*)b)->value;
}

long* reduce_bw_of_points(
	edge_neighbour* edges, 
	long nedges, 
	long n_internal_edges, 
	double* volumes, 
	double3* coords,
	long npoints, 
	long* mg_connectivity, 
	long mg_connectivity_size, 
	long* mg_connectivity_level_below, 
	long mg_connectivity_level_below_size, 
	bool reverse) 
{
	int* point_degree = 
		count_point_internal_degrees(
			npoints, 
			edges, 
			0,
			nedges-1);

	long** point_neighbours = 
		identify_pointwise_internal_neighbours(
			point_degree, 
			npoints, 
			edges, 
			0, 
			nedges-1);

	long* permtab = alloc<long>(npoints);
	for (long i=0; i<npoints; i++) {
		permtab[i] = -1;
	}

	bool* visited_points = alloc<bool>(npoints);
	for (long i=0; i<npoints; i++) {
		visited_points[i] = false;
	}

	// Find lowest-degree point and set to first element of reordered list:
	long lowest_degree_point = 0;
	long lowest_degree = point_degree[0];
	for (long i=1; i<npoints; i++) {
		if (point_degree[i] < lowest_degree) {
			lowest_degree_point = i;
			lowest_degree = point_degree[i];
		}
	}
	long* point_queue = alloc<long>(npoints);
	for (long i=0; i<npoints; i++) point_queue[i] = -1;
	bool* queued_points = alloc<bool>(npoints);
	for (long i=0; i<npoints; i++) queued_points[i] = false;
	long perm_write = 0;
	point_queue[0] = lowest_degree_point;
	queued_points[lowest_degree_point] = true;
	long queue_read=0, queue_write=1;

	show_progress_bar();

	// Apply Cuthill-McKee algorithm to find a reordering:
	long search_start = 0;
	while (queue_read < npoints) {
		long p = point_queue[queue_read++];
		if (p==-1) {
			fprintf(stderr, "ERROR: popped p is -1\n");
			printf("%s:%d\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		if (!queued_points[p]) {
			fprintf(stderr, "ERROR: Have popped point %d that is not marked as queued\n", p);
			printf("%s:%d\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		if (visited_points[p])
			continue;
		permtab[p] = perm_write++;
		visited_points[p] = true;

		struct tagged_int neighbours[point_degree[p]];
		for (long n=0; n<point_degree[p]; n++) {
			neighbours[n].value = point_degree[point_neighbours[p][n]];
			neighbours[n].tag = point_neighbours[p][n];
		}
		qsort(neighbours, point_degree[p], sizeof(tagged_int), compare_two_tagged_ints);

		for (long n=0; n<point_degree[p]; n++) {
			long pn = neighbours[n].tag;
			if (pn >= npoints || pn < 0) {
				fprintf(stderr, "ERROR: pn=%d is out of range [0, %d]\n", pn, npoints-1);
				printf("%s:%d\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			if (!queued_points[pn]) {
				if (queue_write >= npoints) {
					fprintf(stderr, "ERROR: attempting to queue %dth neighbour of %d but queue_write=%d is >= npoints=%d (queue_read=%d)\n", n, p, queue_write, npoints, queue_read);
					fprintf(stderr, "queued_points[%d] = %d\n", pn, queued_points[pn]);
					fprintf(stderr, "visited_points[%d] = %d\n", pn, visited_points[pn]);
					printf("%s:%d\n", __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				point_queue[queue_write++] = pn;
				queued_points[pn] = true;
			}
		}

		if (queue_read%5000==0) {
			update_progress_bar((float)queue_read/(float)npoints);
		}

		if (queue_read == queue_write && queue_read < npoints) {
			// Have reached end of queue but have not processed all
			// points so there must be another disconnected mesh. 
			// Find lowest-degree point of that disconnected mesh 
			// and continue:
			lowest_degree = -1;
			lowest_degree_point = -1;
			for (long p=search_start; p<npoints; p++) {
				if (!visited_points[p] && !queued_points[p]) {
					if ((lowest_degree == -1) || (point_degree[p] < lowest_degree)) {
						lowest_degree_point = p;
						lowest_degree = point_degree[p];
						search_start = p+1;
					}
				}
			}

			if (lowest_degree_point == -1) {
				fprintf(stderr, "ERROR: could not find an unprocessed point\n");
				printf("%s:%d\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			
			point_queue[queue_write++] = lowest_degree_point;
			queued_points[lowest_degree_point] = true;
		}
	}
	dealloc<long>(point_queue);
	dealloc<bool>(queued_points);

	update_progress_bar(1.0);
	stop_progress_bar();

	if (reverse) {
		for (long i=0; i<npoints; i++) {
			permtab[i] = (npoints-1)-permtab[i];
		}
	}

	// Then renumber the a's and b's of edges.
	for (long e=0; e<nedges; e++) {
		if (edges[e].a >= 0)
			edges[e].a = permtab[edges[e].a];
		if (edges[e].b >= 0)
			edges[e].b = permtab[edges[e].b];
	}

	permute_array_inplace(volumes, npoints, permtab);
	permute_array_inplace(coords,  npoints, permtab);

	// Update multigrid connectivity mappings:
	if (mg_connectivity != NULL) {
		// Update LHS of the mapping int_(L) -> int_(L+1), which belongs to this level:
		permute_array_inplace(mg_connectivity, mg_connectivity_size, permtab);
	}
	if (mg_connectivity_level_below != NULL) {
		// Update RHS of the mapping int_(L-1) -> int_(L), which belongs to the level below:
		for (long i=0; i<mg_connectivity_level_below_size; i++) {
			mg_connectivity_level_below[i] = permtab[mg_connectivity_level_below[i]];
		}
	}

	// Create inverse of permtab as 'peritab':
	long* peritab = alloc<long>(npoints);
	for (long p=0; p<npoints; p++) {
		peritab[permtab[p]] = p;
	}
	return peritab;
}
