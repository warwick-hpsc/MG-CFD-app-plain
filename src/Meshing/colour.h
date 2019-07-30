#ifndef COLOUR_H
#define COLOUR_H

#include "common.h"

#define UNCOLOURED -1

struct edge_record { int eidx, coloured_neighbours; };

// Count the number of coloured edge neighbours
int count_coloured_neighbourhood(
	int edgeIdx, 
	int edge_degree, 
	const int* neighbours, 
	const int* edge_colours);

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
	int* number_of_colours);

#endif
