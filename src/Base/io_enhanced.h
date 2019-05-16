#ifndef IO_ENHANCED_H
#define IO_ENHANCED_H

#include "common.h"

std::string generate_input_binary_filename_suffix();

std::string generate_output_filename_suffix(int level);

std::string generate_output_filepath(std::string filename, int level);

std::string generate_solution_filepath(std::string filename, int level);

bool file_exists(const char* filepath);

void duplicate_mesh(
    int* nel,
    double** volumes,
    double3** coords,
    int* number_of_edges,
    int* num_internal_edges, 
    int* num_boundary_edges, 
    int* num_wall_edges, 
    int* boundary_edges_start, 
    int* wall_edges_start,
    edge_neighbour** edges,
    int nel_above,
    int** mg_mapping,
    int* mgc);

bool read_grid_from_bin(
    const char* data_file_name, 
    int* nel, 
    double** volumes, 
    int* number_of_edges, 
    int* num_internal_edges, 
    int* num_boundary_edges, 
    int* num_wall_edges, 
    int* internal_edges_start, 
    int* boundary_edges_start, 
    int* wall_edges_start, 
    edge_neighbour** edges, 
    double3** coords, 
    int** mg_connectivity, 
    int* mg_connectivity_size);

bool write_grid_to_bin(
    const char* data_file_name, 
    int nel, 
    const double *restrict volumes, 
    int number_of_edges, 
    int num_internal_edges, 
    int num_boundary_edges, 
    int num_wall_edges, 
    int internal_edges_start, 
    int boundary_edges_start, 
    int wall_edges_start, 
    const edge_neighbour* edges, 
    const double3* coords, 
    const int* mg_connectivity, 
    int mg_connectivity_size);

void read_input_dat(
    const char* file_name, 
    int* problem_size, 
    std::string** layers, 
    std::string** mg_connectivity);

#ifdef PAPI
void read_papi_config(const char* filename, int* nevents, int** events);
#endif

void read_mg_connectivity(const char* file_name, int** mg_connectivity, int* mgc);

bool read_double_array(
    double* data_array, std::string name, 
    int ndim, int nel, int level);

void dump_step_factors(
    const double *restrict step_factors, 
    int npoints, int level);

void dump_edge_fluxes(
    const edge* edge_variables,
    int num_internal_edges, 
    int num_boundary_edges, 
    int num_wall_edges, 
    int internal_edges_start, 
    int boundary_edges_start, 
    int wall_edges_start, 
    int level);

void dump_flux(
    const double *restrict fluxes, 
    int npoints, int level);

void dump_volumes(
    const double *restrict volumes, 
    int npoints, int level);

void prepare_csv_identification(
    bool write_header,
    std::string *header_s,
    std::string *data_line_s,
    int input_size);

#endif

