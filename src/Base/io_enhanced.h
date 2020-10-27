#ifndef IO_ENHANCED_H
#define IO_ENHANCED_H

#include "common.h"

std::string generate_input_binary_filename_suffix();

std::string generate_output_filename_suffix(int level);

std::string generate_output_filepath(std::string filename, int level);

std::string generate_solution_filepath(std::string filename, int level);

bool file_exists(const char* filepath);

void duplicate_mesh(
    long* nel,
    double** volumes,
    double3** coords,
    long* number_of_edges,
    long* num_internal_edges, 
    long* num_boundary_edges, 
    long* num_wall_edges, 
    long* boundary_edges_start, 
    long* wall_edges_start,
    #if defined BIN_COLOURED_VECTORS || defined BIN_COLOURED_CONTIGUOUS
        long* internal_serial_section_start,
        long* boundary_serial_section_start, 
        long* wall_serial_section_start,
    #endif
    edge_neighbour** edges,
    long nel_above,
    long** mg_mapping,
    long* mgc,
    long** peritab);

bool read_grid_from_bin(
    const char* data_file_name, 
    long* nel, 
    double** volumes, 
    long* number_of_edges, 
    long* num_internal_edges, 
    long* num_boundary_edges, 
    long* num_wall_edges, 
    long* internal_edges_start, 
    long* boundary_edges_start, 
    long* wall_edges_start, 
    edge_neighbour** edges, 
    double3** coords, 
    long** mg_connectivity, 
    long* mg_connectivity_size);

bool write_grid_to_bin(
    const char* data_file_name, 
    long nel, 
    const double *restrict volumes, 
    long number_of_edges, 
    long num_internal_edges, 
    long num_boundary_edges, 
    long num_wall_edges, 
    long internal_edges_start, 
    long boundary_edges_start, 
    long wall_edges_start, 
    const edge_neighbour* edges, 
    const double3* coords, 
    const long* mg_connectivity, 
    long mg_connectivity_size);

void read_input_dat(
    const char* file_name, 
    int* problem_size, 
    std::string** layers, 
    std::string** mg_connectivity);

#ifdef PAPI
void read_papi_config(const char* filename, int* nevents, int** events);
#endif

void read_mg_connectivity(const char* file_name, long** mg_connectivity, long* mgc);

bool read_double_array(
    double* data_array, std::string name, 
    int ndim, long nel, int level);

void dump_step_factors(
    const double *restrict step_factors, 
    long npoints, int level);

void dump_edge_fluxes(
    const edge* edge_variables,
    long num_internal_edges, 
    long num_boundary_edges, 
    long num_wall_edges, 
    long internal_edges_start, 
    long boundary_edges_start, 
    long wall_edges_start, 
    int level);

void dump_flux(
    const double *restrict fluxes, 
    long npoints, int level);

void dump_volumes(
    const double *restrict volumes, 
    long npoints, int level);

void prepare_csv_identification(
    bool write_header,
    std::string *header_s,
    std::string *data_line_s,
    int input_size);

#endif

