#ifndef REDUCE_BW_H
#define REDUCE_BW_H

#include "common.h"

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
	bool reverse);

// Reorder the elements of an array in-place without using any extra
// memory. 'mapping' defines the indices of where elements should be
// moved to, eg mapping[0]=5 means array[0] should be moved to array[5].
template <typename T>
void permute_array_inplace(T* array, long length, long* mapping)
{
    // Reorder in-place
    T temp1, temp2;
    long temp1_idx, temp2_idx;
    long reordered_count=0;
    bool* done = alloc<bool>(length);
    long search_start=0;
    for (long e=0; e<length; e++) {
        if (mapping[e] == e) {
            done[e] = true;
            reordered_count++;
        } else {
            done[e] = false;
        }
    }
    while (reordered_count != length) {
        temp1_idx = -1;
        for (long i=search_start; i<length; i++) {
            if (!done[i]) {
                temp1 = array[i];
                temp1_idx = i;
                search_start=i+1;
                break;
            }
        }
      #ifdef TESTS
        if (temp1_idx == -1) {
            printf("ERROR: Have not reordered all array items but cannot find an unordered item.\n");
            printf("%s:%d\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
      #endif

        while(!done[temp1_idx]) {
          #ifdef TESTS
            if (temp1_idx >= length || temp1_idx < 0) {
                printf("ERROR: temp1_idx=%d is out-of-bounds\n", temp1_idx);
                printf("%s:%d\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            else if (reordered_count > length) {
                printf("ERROR: Have remapped too many elements\n");
                printf("%s:%d\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
          #endif
            // Put 'temp1' into new location as specified by 'mapping',
            // but first make a backup of the destination into 'temp2'.
            long idx_to = mapping[temp1_idx];
          #ifdef TESTS
            if (idx_to == temp1_idx) {
                printf("ERROR: Mapping an item to itself, will happen infinitely\n");
                printf("%s:%d\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            if (idx_to >= length) {
                printf("ERROR: mapping has returned an out-of-range index\n");
                printf("%s:%d\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
          #endif

            temp2 = array[idx_to];
            temp2_idx = idx_to;

          #ifdef TESTS
            if (idx_to >= length) {
                printf("ERROR: mapping[%d] = %d, which is out of bounds.\n", temp1_idx, mapping[temp1_idx]);
                printf("%s:%d\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
          #endif

            array[mapping[temp1_idx]] = temp1;
            done[temp1_idx] = true;
            reordered_count++;

            temp1 = temp2;
            temp1_idx = temp2_idx;
        }
    }
    dealloc<bool>(done);
}

#endif
