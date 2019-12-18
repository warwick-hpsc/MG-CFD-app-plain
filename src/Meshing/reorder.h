#ifndef REORDER_H
#define REORDER_H

#include "common.h"

#ifdef BIN_COLOURED_VECTORS
void BinEdgesIntoColouredVectorUnits(
	long estart, 
	long eend,
	edge_neighbour* edges, 
	long nedges, 
	int* edge_colours,
	int ncolours,
	long nel, 
	long* serial_section_start);
#endif

#ifdef BIN_COLOURED_CONTIGUOUS
void BinEdgesIntoContiguousColouredBlocks(
	long estart, 
	long eend,
	edge_neighbour* edges, 
	long nedges, 
	int* edge_colours,
	int ncolours, 
	long nel,
	long* serial_section_start);
#endif

#endif