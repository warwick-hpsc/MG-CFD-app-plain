// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

#include <iostream>
#include <stdlib.h>

#ifdef LEGACY_ORDERING
#include <algorithm>
#endif

#include "io.h"
#include "io_enhanced.h"

void read_grid(
    const char* data_file_name, 
    long* nel, 
    double** volumes, 
    long* number_of_edges, 
    #ifdef E_SOA
        long* edges_soa_step, 
    #endif
    long* num_internal_edges, 
    long* num_boundary_edges, 
    long* num_wall_edges, 
    long* internal_edges_start, 
    long* boundary_edges_start, 
    long* wall_edges_start, 
    edge_neighbour** edges, 
    double3** coords)
{
    log("read_grid() called");

    // read in domain geometry
    // create edge list
    *number_of_edges = 0;
    *num_internal_edges = 0;
    *num_boundary_edges = 0;
    *num_wall_edges = 0;
    
    log("Opening '%s'\n", data_file_name);

    std::ifstream file(data_file_name);
    if(!file.is_open())
    {
        std::cout << "could not open data file: '" << data_file_name << "'" << std::endl;
        DEBUGGABLE_ABORT
    }

    std::ifstream coords_file((std::string(data_file_name) + std::string(".coords")).c_str());
    if(!coords_file.is_open() && levels > 1)
    {
        std::cout << "could not open coords file for: " << data_file_name << std::endl;
        DEBUGGABLE_ABORT
    }

    file >> (*nel);
    file >> (*number_of_edges);
    edge_neighbour* edge_bin = alloc<edge_neighbour>(*number_of_edges);

    long* internal_edge_indices = alloc<long>(*number_of_edges);
    long* boundary_edge_indices = alloc<long>(*number_of_edges);
    long* wall_edge_indices = alloc<long>(*number_of_edges);

    *volumes = alloc<double>(*nel);
    *coords = alloc<double3>(*nel);
    long** point_neighbours = alloc<long*>(*nel);

    // read in data
    long edge_count=0;
    for(long i = 0; i < *nel; i++)
    {
        file >> (*volumes)[i];
        int point_degree;
        file >> point_degree;
        point_neighbours[i] = alloc<long>(point_degree);

        if (levels > 1) {
            coords_file >> (*coords)[i].x;
            coords_file >> (*coords)[i].y;
            coords_file >> (*coords)[i].z;
        }

        for (int j=0; j < point_degree; j++)
        {
            file >> point_neighbours[i][j];

            double3 edge_weight;
            file >> edge_weight.x;
            file >> edge_weight.y;
            file >> edge_weight.z;

            const long i2 = point_neighbours[i][j];
            if (i2 < i)
            {
                if (i2  == -1)
                {
                    boundary_edge_indices[*num_boundary_edges] = edge_count;
                    (*num_boundary_edges)++;
                }
                else if (i2 == -2)
                {
                    wall_edge_indices[*num_wall_edges] = edge_count;
                    (*num_wall_edges)++;
                }
                else
                {
                    internal_edge_indices[*num_internal_edges] =  edge_count;
                    (*num_internal_edges)++;
                }

                edge_bin[edge_count].a = i2;
                edge_bin[edge_count].b = i;
                edge_bin[edge_count].x = edge_weight.x;
                edge_bin[edge_count].y = edge_weight.y;
                edge_bin[edge_count].z = edge_weight.z;

                if (mesh_variant == MESH_FVCORR) {
                    // rodinia/cfd is flipping all normals, so repeat here
                    edge_bin[edge_count].x *= -1;
                    edge_bin[edge_count].y *= -1;
                    edge_bin[edge_count].z *= -1;
                }
                else {
                    // ... but blindly flipping all normals makes no sense, 
                    // so handle other meshes differently.
                    if (i2 >= 0) {
                        // Is an internal edge, added backwards, so 
                        // need to flip the weight:
                        edge_bin[edge_count].x *= -1;
                        edge_bin[edge_count].y *= -1;
                        edge_bin[edge_count].z *= -1;
                    }
                }

                edge_count++;
            }
        }
    }
    for(long i = 0; i < *nel; i++)
    {
        dealloc<long>(point_neighbours[i]);
    }
    dealloc<long*>(point_neighbours);

    if (edge_count != *number_of_edges) {
        printf("WARNING: Mesh claims to have %ld edges, actually has %ld\n", *number_of_edges, edge_count);
    }

    *internal_edges_start = 0;
    *boundary_edges_start = *num_internal_edges;
    *wall_edges_start     = *boundary_edges_start + *num_boundary_edges;

    *edges = alloc<edge_neighbour>(*number_of_edges);
    long j=0;
    for(long i=0; i < *num_internal_edges; i++,j++)
    {
        const long idx = internal_edge_indices[i];
        (*edges)[j] = edge_bin[idx];
    }
    for (; j<*boundary_edges_start; j++) {
        (*edges)[j].a = -5;
        (*edges)[j].b = -5;
    }
    for(long i=0; i < *num_boundary_edges; i++,j++)
    {
        const long idx = boundary_edge_indices[i];
        (*edges)[j] = edge_bin[idx];
    }
    for (; j<*wall_edges_start; j++) {
        (*edges)[j].a = -5;
        (*edges)[j].b = -5;
    }
    for(long i=0; i < *num_wall_edges; i++,j++)
    {
        const long idx = wall_edge_indices[i];
        (*edges)[j] = edge_bin[idx];
    }   
    for (; j<*number_of_edges; j++) {
        (*edges)[j].a = -5;
        (*edges)[j].b = -5;
    }

    #ifdef LEGACY_ORDERING
        std::sort(*edges, 
                  (*edges)+(*num_internal_edges), 
                  compare_two_edges);
        std::sort(*edges+(*num_internal_edges), 
                  (*edges)+(*num_internal_edges)+(*num_boundary_edges), 
                  compare_two_edges);
        std::sort(*edges+(*num_internal_edges)+(*num_boundary_edges), 
                  (*edges)+(*num_internal_edges)+(*num_boundary_edges)+(*num_wall_edges), 
                  compare_two_edges);
    #endif

    dealloc<edge_neighbour>(edge_bin);
    dealloc<long>(internal_edge_indices);
    dealloc<long>(boundary_edge_indices);
    dealloc<long>(wall_edge_indices);
}

void dump(
    const double *restrict variables, 
    long nel, 
    int level)
{
    log("dump_variables() called");

    std::string filepath = generate_output_filepath(std::string("variables"), level);
    FILE *file = fopen(filepath.c_str(), "w");
    if (file == NULL) {
        fprintf(stderr, "ERROR: Failed to open file for writing: '%s'\n", filepath.c_str());
        exit(EXIT_FAILURE);
    }

    printf("Dumping variables[] to file: %s\n", filepath.c_str());

    for (long i=0; i<nel; i++) {
        const long p_idx  = i*NVAR + VAR_DENSITY;
        const long mx_idx = i*NVAR + VAR_MOMENTUMX;
        const long my_idx = i*NVAR + VAR_MOMENTUMY;
        const long mz_idx = i*NVAR + VAR_MOMENTUMZ;
        const long pe_idx = i*NVAR + VAR_DENSITY_ENERGY;

        fprintf(file, "%.17e %.17e %.17e %.17e %.17e\n", variables[p_idx], 
                                                         variables[mx_idx], 
                                                         variables[my_idx], 
                                                         variables[mz_idx],
                                                         variables[pe_idx]);
    }
    fclose(file);
    
    log("dump() complete");
}
