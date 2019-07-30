#ifndef REORDER_H
#define REORDER_H

#include "common.h"

#ifdef BIN_COLOURED_VECTORS
void BinEdgesIntoColouredVectorUnits(
	int estart, 
	int eend,
	edge_neighbour* edges, 
	int nedges, 
	int* edge_colours,
	int ncolours,
	int nel, 
	int* serial_section_start);
#endif

#ifdef BIN_COLOURED_CONTIGUOUS
void BinEdgesIntoContiguousColouredBlocks(
	int estart, 
	int eend,
	edge_neighbour* edges, 
	int nedges, 
	int* edge_colours,
	int ncolours, 
	int nel,
	int* serial_section_start);
#endif

#endif