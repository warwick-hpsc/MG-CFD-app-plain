#include "reorder.h"

bool validate_permutation(
	long* permtab,
	long map_size)
{
	// Test that the proposed permuation is valid:
	bool success=true;
	bool* post_permute_record = alloc<bool>(map_size);
	for (long i=0; i<map_size; i++) post_permute_record[i] = false;

	for (long i=0; i<map_size; i++) {
		long to_idx = permtab[i];
		if (to_idx >= map_size) {
			fprintf(stderr, "ERROR: permtab[%d] -> to_idx=%d is >= map_size=%d\n", i, to_idx, map_size);
			success = false;
			break;
		}
		if (to_idx < 0) {
			fprintf(stderr, "ERROR: permtab[%d] -> to_idx=%d is < 0\n", i, to_idx);
			success = false;
			break;
		}
		if (post_permute_record[to_idx]) {
			fprintf(stderr, "ERROR: Permtab maps %d to %d, but a previous item has already been mapped to that location\n", i, to_idx);
			success = false;
			break;
		}
		post_permute_record[to_idx] = true;
	}
	if (success) {
		for (long i=0; i<map_size; i++) {
			if (!post_permute_record[i]) {
				fprintf(stderr, "ERROR: No item has been mapped to %d\n", i);
				success = false;
				break;
			}
		}
	}
	dealloc<bool>(post_permute_record);

	if (!success) {
		fprintf(stderr, "mapsize = %d\n", map_size);
	}

	return success;
}

#ifdef BIN_COLOURED_VECTORS
void BinEdgesIntoColouredVectorUnits(
	long estart, 
	long eend,
	edge_neighbour* edges, 
	long nedges, 
	int* edge_colours,
	int ncolours, 
	long nel,
	long* serial_section_start)
{
	log("BinEdgesIntoColouredVectorUnits()");

	#ifdef TESTS
		if (eend < estart) {
			fprintf(stderr, "ERROR: eend=%d is < estart=%d\n", eend, estart);
			DEBUGGABLE_ABORT
		}
	#endif

	if (DBLS_PER_SIMD == 1) {
		// Hardcode this specific case, primarily because there is a 
		// bug in the logic below triggered only when DBLS_PER_SIMD
		// is 1.
		*serial_section_start = eend+1;
		return;
	}

	long total_range = eend-estart+1;
	long* permtab = alloc<long>(total_range);
	for (long i=0; i<total_range; i++) permtab[i] = -1;
	long next_permtab_idx = 0;

	long num_processed_edges = 0;
	bool* processed_edges = alloc<bool>(total_range);
	for (long i=0; i<total_range; i++) processed_edges[i] = false;
	long next_edge = estart;
	long vector_block[DBLS_PER_SIMD];
	int vector_colour;
	int vector_idx = 0;
	while ((next_edge-estart) <= total_range) {
		// Find an edge that has not been processed:
		while(processed_edges[next_edge-estart]) {
			next_edge++;
		}

		// Specify colour of this vector block:
		vector_colour = edge_colours[next_edge];

		// Fill vector block:
		vector_idx = 0;
		vector_block[vector_idx] = next_edge;
		vector_idx++;
		next_edge++;
		for (long e=next_edge; (vector_idx<DBLS_PER_SIMD) && (e<=eend); e++) {
			if (edge_colours[e] == vector_colour) {
				vector_block[vector_idx] = e;
				vector_idx++;
			}
		}

		if (vector_idx == DBLS_PER_SIMD) {
			// Construct a permutation for this vector block, so that
			// after application the edges are contiguous:
			for (int i=0; i<vector_idx; i++) {
				long e = vector_block[i];
				permtab[e-estart] = next_permtab_idx;
				next_permtab_idx++;

				processed_edges[e-estart] = true;
				num_processed_edges++;
			}
		}
	}
	if (num_processed_edges != total_range) {
		*serial_section_start = estart + next_permtab_idx;
		for (long e=estart; e<=eend; e++) {
			if (!processed_edges[e-estart]) {
				permtab[e-estart] = next_permtab_idx;
				next_permtab_idx++;
			}
		}
	} else {
		// No serial section.
		*serial_section_start = eend+1;
	}

    if (!validate_permutation(permtab, total_range)) {
    	DEBUGGABLE_ABORT
    }

	// Apply the reordering:
	permute_array_range(edges,        permtab, estart, total_range);
	permute_array_range(edge_colours, permtab, estart, total_range);

	dealloc<bool>(processed_edges);
	processed_edges = NULL;

	log("BinEdgesIntoColouredVectorUnits() complete");
}
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
	long* serial_section_start)
{
	log("BinEdgesIntoContiguousColouredBlocks()");

	long total_range = eend-estart+1;
	long* permtab = alloc<long>(total_range);
	for (long i=0; i<total_range; i++) permtab[i] = -1;
	long next_permtab_idx = 0;

	long num_processed_edges = 0;
	bool* processed_edges = alloc<bool>(total_range);
	for (long i=0; i<total_range; i++) processed_edges[i] = false;

	long next_edge = estart;
	long vector_block[DBLS_PER_SIMD];
	int vector_colour;
	int vector_idx = 0;
	for (int c=0; c<ncolours; c++) {
		// Find an edge that has not been processed:
		next_edge = estart;
		while(next_edge <= eend && edge_colours[next_edge] != c) {
			next_edge++;
		}
		if (next_edge > eend) {
			// Did not find any edges with colour 'c'
			continue;
		}

		while(true) {
			vector_block[0] = next_edge;
			vector_idx = 1;
			for (long e=next_edge+1; (vector_idx<DBLS_PER_SIMD) && (e<=eend); e++) {
				if (edge_colours[e] == c) {
					vector_block[vector_idx] = e;
					vector_idx++;
				}
			}

			if (vector_idx < DBLS_PER_SIMD) {
				// No more vector-sized edge blocks of colour 'c', 
				// so move onto next colour:
				break;
			}

			if (vector_idx == DBLS_PER_SIMD) {
				// Construct a permutation for this vector block, so that
				// after application the edges are contiguous:
				for (int i=0; i<vector_idx; i++) {
					long e = vector_block[i];
					// permtab[e-estart] = estart + next_permtab_idx;
					permtab[e-estart] = next_permtab_idx;
					next_permtab_idx++;

					processed_edges[e-estart] = true;
					num_processed_edges++;
				}
			}

			next_edge = vector_block[DBLS_PER_SIMD-1] + 1;
			while(next_edge <= eend && edge_colours[next_edge] != c) {
				next_edge++;
			}
			if (next_edge > eend) {
				break;
			}
		}
	}

	if (num_processed_edges != total_range) {
		*serial_section_start = estart + next_permtab_idx;
		for (long e=estart; e<=eend; e++) {
			if (!processed_edges[e-estart]) {
				permtab[e-estart] = next_permtab_idx;
				next_permtab_idx++;
			}
		}
	} else {
		// No serial section.
		*serial_section_start = eend+1;
	}

	// #ifdef TESTS
	    if (!validate_permutation(permtab, total_range)) {
	    	DEBUGGABLE_ABORT
	    }
    // #endif

	// Apply the reordering:
	permute_array_range(edges,        permtab, estart, total_range);
	permute_array_range(edge_colours, permtab, estart, total_range);

	dealloc<bool>(processed_edges);
	processed_edges = NULL;

	log("BinEdgesIntoContiguousColouredBlocks() complete");
}
#endif
