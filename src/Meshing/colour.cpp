#include <queue>
#include <iterator>
#include <stdlib.h>
// #include <string.h>

#include "colour.h"
#include "graph.h"
#include "progress.h"

class CompareEdgeRecords {
public:
    bool operator()(const edge_record& e1, const edge_record& e2) {
        return e1.coloured_neighbours < e2.coloured_neighbours;
    }
};

// Count the number of coloured edge neighbours
int count_coloured_neighbourhood(
	int edgeIdx, 
	int edge_degree, 
	const int* neighbours, 
	const int* edge_colours)
{
	int count = 0;
	for (int i=0; i<edge_degree; i++) {
		if (edge_colours[neighbours[i]] != UNCOLOURED) {
			count++;
		}
	}
	return count;
}

// Colour the edges of mesh so that no two connected edges have the same colour.
void colour_mesh_strict(
	const edge_neighbour* edges, 
	int nedges, 
	int npoints, 
    int boundary_edges_start, 
    int wall_edges_start, 
    int num_internal_edges, 
    int num_boundary_edges, 
    int num_wall_edges, 
	int** edge_colours, 
	int* number_of_colours)
{
	// Prepare arrays that contain mesh connectivity information:
	int* point_degree = 
		count_point_degrees(
			npoints, 
			edges, 
			0, 
			nedges-1);

	int** point_to_edges = 
		construct_point_to_edge_mapping(
			point_degree, 
			npoints, 
			edges, 
			0, 
			nedges-1);

	int* edge_degree = 
		count_edge_degrees(
			point_degree, 
			point_to_edges, 
			edges, 
			nedges, 
			0, 
			nedges-1);

	int** edge_neighbours =
		identify_edgewise_neighbours(
			point_to_edges, 
			point_degree, 
			edge_degree, 
			edges, 
			nedges, 
			0, 
			nedges-1);

	int max_edge_degree = 0;
	for (int e=0; e<nedges; e++) {
		int d = edge_degree[e];
		if (d > max_edge_degree) max_edge_degree = d;
	}

	show_progress_bar();

	// Colour the edges of mesh so that no two connected edges have the same colour.
	*edge_colours = alloc<int>(nedges);
	bool* colour_availability = alloc<bool>(max_edge_degree);
	for (int i=0; i<nedges; i++) (*edge_colours)[i] = UNCOLOURED;
	bool uncoloured_edges_remain=true;

	int* coloured_neighbour_counts = alloc<int>(nedges);
	for (int i=0; i<nedges; i++) coloured_neighbour_counts[i] = 0;

	std::priority_queue<edge_record, std::vector<edge_record>, CompareEdgeRecords> edge_queue;
	edge_record first;
	first.eidx=0;
	first.coloured_neighbours=0;
	edge_queue.push(first);
	int nedges_coloured=0;
	int search_start=0;
	while(uncoloured_edges_remain) {
		while (!edge_queue.empty()) {
			int eid = edge_queue.top().eidx;
			edge_queue.pop();
			if ((*edge_colours)[eid] != UNCOLOURED) {
				// Edge has already been coloured
				continue;
			}

			for (int i=0; i<max_edge_degree; i++) colour_availability[i] = true;

			for (int n=0; n<edge_degree[eid]; n++) {
				int neighbour_colour = (*edge_colours)[edge_neighbours[eid][n]];
				if (neighbour_colour != UNCOLOURED) 
					colour_availability[neighbour_colour] = false;
			}
			for (int c=0; c<max_edge_degree; c++) {
				if (colour_availability[c] == true) {
					(*edge_colours)[eid] = c;
					break;
				}
			}
			nedges_coloured++;
			if (nedges_coloured%5000==0) {
				update_progress_bar((float)nedges_coloured/(float)nedges);
			}
			// Push on edge-neighbours of newly-coloured edge.
			for (int n=0; n<edge_degree[eid]; n++) {
				coloured_neighbour_counts[edge_neighbours[eid][n]]++;
				edge_record er;
				er.eidx = edge_neighbours[eid][n];
				er.coloured_neighbours = coloured_neighbour_counts[edge_neighbours[eid][n]];
				edge_queue.push(er);
			}
		}

		uncoloured_edges_remain = false;

		// Search for next uncoloured edge, skipping any padding edges:
		for (int e=search_start; e<num_internal_edges; e++) {
			if ((*edge_colours)[e] == UNCOLOURED) {
				uncoloured_edges_remain = true;
				edge_record er;
				er.eidx = e;
				er.coloured_neighbours = coloured_neighbour_counts[e];
				edge_queue.push(er);
				search_start=e+1;
				break;
			}
		}
		if (!uncoloured_edges_remain) {
			search_start = boundary_edges_start;
			for (int e=search_start; e<boundary_edges_start+num_boundary_edges; e++) {
				if ((*edge_colours)[e] == UNCOLOURED) {
					uncoloured_edges_remain = true;
					edge_record er;
					er.eidx = e;
					er.coloured_neighbours = coloured_neighbour_counts[e];
					edge_queue.push(er);
					search_start=e+1;
					break;
				}
			}
		}
		if (!uncoloured_edges_remain) {
			search_start = wall_edges_start;
			for (int e=search_start; e<wall_edges_start+num_wall_edges; e++) {
				if ((*edge_colours)[e] == UNCOLOURED) {
					uncoloured_edges_remain = true;
					edge_record er;
					er.eidx = e;
					er.coloured_neighbours = coloured_neighbour_counts[e];
					edge_queue.push(er);
					search_start=e+1;
					break;
				}
			}
		}
	}

	update_progress_bar(1.0);
	stop_progress_bar();

	// Check all edges have been coloured
    #ifdef TESTS
		for (int e=0; e<num_internal_edges; e++) {
			if ((*edge_colours)[e] == UNCOLOURED) {
				printf("Edge %d not coloured. a=%d, b=%d\n", e, edges[e].a, edges[e].b);
				// Does edge have neighbours?
				for (int en=0; en<nedges; en++) {
					int a=edges[e].a, b=edges[e].b, an=edges[en].a, bn=edges[en].b;
					if ( (a>=0 && ((a==an && b!=bn) || (a==bn && b!=an))) ||
						 (b>=0 && ((b==an && a!=bn) || (b==bn && a!=an))) ) {
						printf("Edge %d is a neighbour of %d\n", en, e);
					}
				}
			}
		}
		for (int e=boundary_edges_start; e<boundary_edges_start+num_boundary_edges; e++) {
			if ((*edge_colours)[e] == UNCOLOURED) {
				printf("Edge %d not coloured. a=%d, b=%d\n", e, edges[e].a, edges[e].b);
				// Does edge have neighbours?
				for (int en=0; en<nedges; en++) {
					int a=edges[e].a, b=edges[e].b, an=edges[en].a, bn=edges[en].b;
					if ( (a>=0 && ((a==an && b!=bn) || (a==bn && b!=an))) ||
						 (b>=0 && ((b==an && a!=bn) || (b==bn && a!=an))) ) {
						printf("Edge %d is a neighbour of %d\n", en, e);
					}
				}
			}
		}
		for (int e=wall_edges_start; e<wall_edges_start+num_wall_edges; e++) {
			if ((*edge_colours)[e] == UNCOLOURED) {
				printf("Edge %d not coloured. a=%d, b=%d\n", e, edges[e].a, edges[e].b);
				// Does edge have neighbours?
				for (int en=0; en<nedges; en++) {
					int a=edges[e].a, b=edges[e].b, an=edges[en].a, bn=edges[en].b;
					if ( (a>=0 && ((a==an && b!=bn) || (a==bn && b!=an))) ||
						 (b>=0 && ((b==an && a!=bn) || (b==bn && a!=an))) ) {
						printf("Edge %d is a neighbour of %d\n", en, e);
					}
				}
			}
		}
    #endif

	#ifdef TESTS
		// Validate edge_colours:
		printf("Validating 'edge_colours':\n");
		bool fail=false;
		int* c_counts = alloc<int>(max_edge_degree);
		for (int i=0; i<max_edge_degree; i++) c_counts[i] = 0;
		for (int e=0; e<nedges; e++) {
			if ((*edge_colours)[e] < 0 || (*edge_colours)[e] >= max_edge_degree) {
				printf("Colour of edge %d is out-of-range = %d\n", e, (*edge_colours)[e]);
				fail=true;
			}
			c_counts[(*edge_colours)[e]]++;
		}
		if (!fail)
			printf("Success\n");
		else {
			printf("%s:%d\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		printf("\n");
		dealloc<int>(c_counts);
	#endif

	*number_of_colours = max_edge_degree;

	for (int p=0; p<npoints; p++) {
		if (point_degree[p] > 0)
			dealloc<int>(point_to_edges[p]);
	}
	dealloc<int*>(point_to_edges);
	dealloc<int>(point_degree);
	for (int e=0; e<nedges; e++) {
		if (edge_degree[e] > 0)
			dealloc<int>(edge_neighbours[e]);
	}
	dealloc<int*>(edge_neighbours);
	dealloc<int>(edge_degree);
	dealloc<bool>(colour_availability);
	dealloc<int>(coloured_neighbour_counts);
}
