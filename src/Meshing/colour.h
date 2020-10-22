#ifndef COLOUR_H
#define COLOUR_H

#include "common.h"

#define UNCOLOURED -1

struct edge_record { long eidx, coloured_neighbours; };

// Colour the edges of mesh so that no two connected edges have the same colour.
void colour_mesh_strict(
	const edge_neighbour* edges, 
	long nedges, 
	long npoints, 
	long boundary_edges_start, 
	long wall_edges_start, 
	long num_internal_edges, 
	long num_boundary_edges, 
	long num_wall_edges, 
	int** edge_colours, 
	int* number_of_colours);

#endif
