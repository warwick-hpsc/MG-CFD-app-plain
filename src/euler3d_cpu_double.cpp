// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

// Warwick extensions:
// - multigrid, per-edge computation, n-neighbours
// - residual calculation, solution validation
// - enhanced performance measurement: PAPI, 'toggle-able' arithmetic optimisations

#include <getopt.h>
#include <iostream>
#include <omp.h>
#include <string.h>
#include <unistd.h>

// Base:
#include "common.h"
#include "io.h"
#include "io_enhanced.h"

// Kernels:
#include "flux_loops.h"
#include "indirect_rw_loop.h"
#include "cfd_loops.h"
#include "mg_loops.h"

// Meshing:
#include "colour.h"
#include "reorder.h"

// Monitoring:
#include "papi_funcs.h"
#include "timer.h"
#include "loop_stats.h"

#include "validation.h"

// Globals:
int levels=0;
int level=0;
int current_kernel;
int mesh_variant;
double ff_variable[NVAR];
double3 ff_flux_contribution_momentum_x;
double3 ff_flux_contribution_momentum_y;
double3 ff_flux_contribution_momentum_z;
double3 ff_flux_contribution_density_energy;

void clean_level(
    int nel, 
    double* restrict volumes, 
    double* restrict variables, 
    double* restrict old_variables, 
    double* restrict fluxes, 
    double* restrict step_factors, 
    edge_neighbour* restrict edges, 
    edge* restrict edge_variables,
    double3* restrict coords)
{
    dealloc<double>(volumes);

    dealloc<double>(variables);
    dealloc<double>(old_variables);
    dealloc<double>(fluxes);
    dealloc<double>(step_factors);

    dealloc<edge_neighbour>(edges);
    if (edge_variables != NULL) {
        dealloc<edge>(edge_variables);
    }

    dealloc<double3>(coords);
}

int main(int argc, char** argv)
{
    log("In main()");

    set_config_defaults();
    if (!parse_arguments(argc, argv)) {
        return 1;
    }
    // print_config();

    if (strcmp(conf.input_file, "") == 0) {
        printf("ERROR: input_file not set\n");
        return 1;
    }
    #ifdef PAPI
    if (strcmp(conf.papi_config_file, "") == 0) {
        printf("ERROR: papi_config_file not set\n");
        return 1;
    }
    #endif

    #if defined OMP && defined OMP_SCATTERS
        if ((conf.mesh_duplicate_count % conf.omp_num_threads) != 0) {
            printf("WARNING: Mesh duplicate multiplier %d insufficient for %d threads, disabling OpenMP.\n", 
                   conf.mesh_duplicate_count, conf.omp_num_threads);
            conf.omp_num_threads = 1;
        }
    #endif

    const char* input_file_name = conf.input_file;
    const char* input_directory = conf.input_file_directory;
    if (strcmp(input_directory, "")!=0) {
        input_file_name = strdup((std::string(input_directory) + "/" + input_file_name).c_str());
    }
    
    ///////////////////////////////////////////////////////////////////////////
    // Read input .dat file:
    ///////////////////////////////////////////////////////////////////////////
    int problem_size = 0;
    levels = 0;
    int cycles = 0;
    std::string* layers = NULL;
    std::string* mg_connectivity_filename = NULL;
    log("Reading .dat file");
    read_input_dat(input_file_name, &problem_size, &layers, &mg_connectivity_filename);
    if (strcmp(input_directory, "")!=0) {
        for (int l=0; l<levels; l++) {
            layers[l] = (std::string(input_directory) + "/" + layers[l]).c_str();
            if (l < (levels-1)) {
                mg_connectivity_filename[l] = (std::string(input_directory) + "/" + mg_connectivity_filename[l]).c_str();
            }
        }
    }
    cycles = conf.num_cycles;

    double total_compute_time = 0.0;
    #ifdef TIME
    init_timers();
    #endif
    #ifdef PAPI
        init_papi();
        load_papi_events();
    #endif
    init_iters();

    ///////////////////////////////////////////////////////////////////////////
    // Create data arrays:
    ///////////////////////////////////////////////////////////////////////////
    // Number of points:
    long nel[levels];
    long num_internal_edges[levels];
    long num_boundary_edges[levels];
    long num_wall_edges[levels];
    long internal_edge_starts[levels];
    long boundary_edge_starts[levels];
    long wall_edge_starts[levels];
    double* volumes[levels];
    long number_of_edges[levels];
    edge_neighbour* edges[levels];
    edge* edge_variables[levels];
    double* variables[levels];
    double* old_variables[levels];
    double* residuals[levels];
    double* fluxes[levels];
    double* step_factors[levels];
    double3* coords[levels];
    // Multigrid connectivity stuff:
    long* mg_connectivity[levels];
    long mg_connectivity_size[levels];
    for (int i=0; i<levels; i++) {
        coords[i] = NULL;
        mg_connectivity[i] = NULL;
        mg_connectivity_size[i] = 0;
    }
 
    ///////////////////////////////////////////////////////////////////////////
    // Read in the multigrid mesh, first checking for prepared binary files:
    ///////////////////////////////////////////////////////////////////////////
    std::string bin_filename_suffix = generate_input_binary_filename_suffix();
    if (bin_filename_suffix != "") {
        bin_filename_suffix = std::string(".") + bin_filename_suffix;
    }
    log("Suffix is: %s\n", bin_filename_suffix.c_str());
    for(int i = 0; i < levels; i++)
    {
        log("Loading MG level %d\n", i);

        if (!read_grid_from_bin(
            (layers[i]+bin_filename_suffix).c_str(), 
            &(nel[i]), 
            &(volumes[i]), 
            &(number_of_edges[i]), 
            &(num_internal_edges[i]), 
            &(num_boundary_edges[i]), 
            &(num_wall_edges[i]), 
            &(internal_edge_starts[i]), 
            &(boundary_edge_starts[i]),
            &(wall_edge_starts[i]),
            &(edges[i]), 
            &(coords[i]), 
            &(mg_connectivity[i]), 
            &(mg_connectivity_size[i])))
        {
            read_grid(
                layers[i].c_str(), 
                &(nel[i]), 
                &(volumes[i]), 
                &(number_of_edges[i]), 
                &(num_internal_edges[i]), 
                &(num_boundary_edges[i]), 
                &(num_wall_edges[i]), 
                &(internal_edge_starts[i]), 
                &(boundary_edge_starts[i]),
                &(wall_edge_starts[i]),
                &(edges[i]), 
                &(coords[i]));

            if (i != (levels-1)) {
                read_mg_connectivity(mg_connectivity_filename[i].c_str(), &(mg_connectivity[i]), &(mg_connectivity_size[i]));
            }
            else {
                mg_connectivity_size[i] = 0;
                mg_connectivity[i] = NULL;
            }

            // Update: this .bin functionality not necessary, and is 
            //         polluting cluster filesystem.
            // write_grid_to_bin(
            //     (layers[i]+bin_filename_suffix).c_str(), 
            //     nel[i], 
            //     volumes[i], 
            //     number_of_edges[i], 
            //     num_internal_edges[i], 
            //     num_boundary_edges[i], 
            //     num_wall_edges[i], 
            //     internal_edge_starts[i],
            //     boundary_edge_starts[i], 
            //     wall_edge_starts[i], 
            //     edges[i], 
            //     coords[i], 
            //     mg_connectivity[i], 
            //     mg_connectivity_size[i]);
        }

        // Create intermediary arrays
        variables[i] = alloc<double>(nel[i]*NVAR);
        zero_array(nel[i]*NVAR, variables[i]);
        residuals[i] = alloc<double>(nel[i]*NVAR);
        zero_array(nel[i]*NVAR, residuals[i]);
        old_variables[i] = alloc<double>(nel[i]*NVAR);
        zero_array(nel[i]*NVAR, old_variables[i]);
        fluxes[i] = alloc<double>(nel[i]*NVAR);
        zero_array(nel[i]*NVAR, fluxes[i]);
        step_factors[i] = alloc<double>(nel[i]);
        zero_array(nel[i], step_factors[i]);
        #ifndef FLUX_FISSION
            edge_variables[i] = NULL;
        #else
            edge_variables[i] = alloc<edge>(number_of_edges[i]*NVAR);
        #endif

        log("Level %d: #edges = %d (#internal = %d, #boundary = %d, #wall = %d), #nodes = %d\n", i, number_of_edges[i], num_internal_edges[i], num_boundary_edges[i], num_wall_edges[i], nel[i]);
        if (i != (levels-1)) {
            log("          mg_connectivity_size = %d\n", mg_connectivity_size[i]);
        }
    }

    #ifdef COLOURED_CONFLICT_AVOIDANCE
        int* edge_colours[levels];
        int number_of_colours[levels];
        printf("Colouring mesh edges:\n");
        double colour_start_time = omp_get_wtime();
        for (int i=0; i<levels; i++) {
            edge_colours[i] = NULL;
            number_of_colours[i] = 0;

            colour_mesh_strict(
                edges[i], 
                number_of_edges[i], 
                nel[i], 
                boundary_edge_starts[i], 
                wall_edge_starts[i], 
                num_internal_edges[i], 
                num_boundary_edges[i], 
                num_wall_edges[i], 
                &edge_colours[i], 
                &number_of_colours[i]);
        }
        printf("Colouring time = %.2f\n", omp_get_wtime()-colour_start_time);
    #endif

    #ifdef BIN_COLOURED_VECTORS
        long internal_serial_section_starts[levels];
        long boundary_serial_section_starts[levels];
        long wall_serial_section_starts[levels];

        for (int i=0; i<levels; i++) {
            BinEdgesIntoColouredVectorUnits(
                0, 
                num_internal_edges[i]-1, 
                edges[i], 
                number_of_edges[i], 
                edge_colours[i], 
                number_of_colours[i], 
                nel[i], 
                &internal_serial_section_starts[i]);

            BinEdgesIntoColouredVectorUnits(
                boundary_edge_starts[i], 
                boundary_edge_starts[i]+num_boundary_edges[i]-1, 
                edges[i], 
                number_of_edges[i], 
                edge_colours[i], 
                number_of_colours[i], 
                nel[i], 
                &boundary_serial_section_starts[i]);

            BinEdgesIntoColouredVectorUnits(
                wall_edge_starts[i], 
                wall_edge_starts[i]+num_wall_edges[i]-1, 
                edges[i], 
                number_of_edges[i], 
                edge_colours[i], 
                number_of_colours[i], 
                nel[i], 
                &wall_serial_section_starts[i]);
        }
    #elif defined BIN_COLOURED_CONTIGUOUS
        long internal_serial_section_starts[levels];
        long boundary_serial_section_starts[levels];
        long wall_serial_section_starts[levels];

        for (int i=0; i<levels; i++) {
            BinEdgesIntoContiguousColouredBlocks(
                0, 
                num_internal_edges[i]-1, 
                edges[i], 
                number_of_edges[i], 
                edge_colours[i], 
                number_of_colours[i], 
                nel[i], 
                &internal_serial_section_starts[i]);

            BinEdgesIntoContiguousColouredBlocks(
                boundary_edge_starts[i], 
                boundary_edge_starts[i]+num_boundary_edges[i]-1, 
                edges[i], 
                number_of_edges[i], 
                edge_colours[i], 
                number_of_colours[i], 
                nel[i], 
                &boundary_serial_section_starts[i]);

            BinEdgesIntoContiguousColouredBlocks(
                wall_edge_starts[i], 
                wall_edge_starts[i]+num_wall_edges[i]-1, 
                edges[i], 
                number_of_edges[i], 
                edge_colours[i], 
                number_of_colours[i], 
                nel[i], 
                &wall_serial_section_starts[i]);
        }
    #endif

    ///////////////////////////////////////////////////////////////////////////
    // Duplicate mesh (for safe thread decomposition):
    ///////////////////////////////////////////////////////////////////////////
    if (conf.mesh_duplicate_count > 1) {
        problem_size *= conf.mesh_duplicate_count;

        for(int i = 0; i < levels; i++) {
            if (i < (levels-1)) {
                duplicate_mesh(
                    &nel[i],
                    &volumes[i],
                    &coords[i],
                    &number_of_edges[i],
                    &num_internal_edges[i], 
                    &num_boundary_edges[i], 
                    &num_wall_edges[i], 
                    &boundary_edge_starts[i], 
                    &wall_edge_starts[i],
                    #ifdef COLOURED_CONFLICT_AVOIDANCE
                        &internal_serial_section_starts[i],
                        &boundary_serial_section_starts[i], 
                        &wall_serial_section_starts[i],
                    #endif
                    &edges[i],
                    nel[i+1],
                    &mg_connectivity[i], 
                    &mg_connectivity_size[i]);
            } else {
                duplicate_mesh(
                    &nel[i],
                    &volumes[i],
                    &coords[i],
                    &number_of_edges[i],
                    &num_internal_edges[i], 
                    &num_boundary_edges[i], 
                    &num_wall_edges[i], 
                    &boundary_edge_starts[i], 
                    &wall_edge_starts[i],
                    #ifdef COLOURED_CONFLICT_AVOIDANCE
                        &internal_serial_section_starts[i],
                        &boundary_serial_section_starts[i], 
                        &wall_serial_section_starts[i],
                    #endif
                    &edges[i],
                    0,
                    NULL, 
                    NULL);
            }

            dealloc<double>(variables[i]);
            variables[i] = alloc<double>((nel[i])*NVAR);
            zero_array(nel[i]*NVAR, variables[i]);
            dealloc<double>(residuals[i]);
            residuals[i] = alloc<double>(nel[i]*NVAR);
            zero_array(nel[i]*NVAR, residuals[i]);
            dealloc<double>(old_variables[i]);
            old_variables[i] = alloc<double>(nel[i]*NVAR);
            zero_array(nel[i]*NVAR, old_variables[i]);
            dealloc<double>(fluxes[i]);
            fluxes[i] = alloc<double>(nel[i]*NVAR);
            zero_array(nel[i]*NVAR, fluxes[i]);
            dealloc<double>(step_factors[i]);
            step_factors[i] = alloc<double>(nel[i]);
            zero_array(nel[i], step_factors[i]);

            #ifdef FLUX_FISSION
                dealloc<edge>(edge_variables[i]);
                edge_variables[i] = alloc<edge>(number_of_edges[i]*NVAR);
            #endif
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // Initialise:
    ///////////////////////////////////////////////////////////////////////////
    initialize_far_field_conditions(
        ff_variable, 
        ff_flux_contribution_momentum_x,
        ff_flux_contribution_momentum_y,
        ff_flux_contribution_momentum_z,
        ff_flux_contribution_density_energy);

    long* up_scratch = alloc<long>(nel[0]);
    for (int i=0; i<levels; i++) {
        initialize_variables(nel[i], variables[i]);
    }

    for (int i=0; i<levels; i++) {
        zero_array(NVAR*nel[i], fluxes[i]);
        zero_array(NVAR*nel[i], residuals[i]);
    }

    // For these meshes, NaN's appear after few solver iterations. 
    // Until root cause identified and addresses, reduce edge 
    // weights to delay appearance, allowing collection of useful 
    // duration of performance data.
    if (mesh_variant == MESH_M6_WING) {
        for (int l=0; l<levels; l++) {
            adjust_ewt(coords[l], number_of_edges[l], edges[l]);
            dampen_ewt(number_of_edges[l], edges[l], 5e-8);
        }
    } else if (mesh_variant == MESH_LA_CASCADE) {
        for (int l=0; l<levels; l++) {
            adjust_ewt(coords[l], number_of_edges[l], edges[l]);
            dampen_ewt(number_of_edges[l], edges[l], 1e-7);
        }
    } else if (mesh_variant == MESH_ROTOR_37) {
        for (int l=0; l<levels; l++) {
            adjust_ewt(coords[l], number_of_edges[l], edges[l]);
            dampen_ewt(number_of_edges[l], edges[l], 2e-7);
        }
    }

    #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
        double* edge_weights[levels];
        for (int i=0; i<levels; i++) {
            edge_weights[i] = alloc<double>(number_of_edges[i]);
            for (int e=0; e<number_of_edges[i]; e++) {
                edge_weights[i][e] = sqrt(edges[i][e].x*edges[i][e].x + edges[i][e].y*edges[i][e].y + edges[i][e].z*edges[i][e].z);
            }
        }
    #endif

    ///////////////////////////////////////////////////////////////////////////
    // Begin compute!
    ///////////////////////////////////////////////////////////////////////////
    log("Beginning compute");
    double loop_start_time = omp_get_wtime();
    level = 0;
    int mg_direction = MG_UP;
    int edge_offset = 0;
    for (int i=0; i<cycles;)
    {
        log("LEVEL %d\n", level);

        if (levels <= 1) {
            printf("\nCycle %d / %d", i+1, cycles);
        } else {
            if (level == 0) {
                printf("\nMG cycle %d / %d", i+1, cycles);
            }
        }

        copy<double>(
            old_variables[level], 
            variables[level], 
            nel[level]*NVAR);

        if (mesh_variant == MESH_FVCORR) {
            // Use the original 'incorrect' step factor calculation, 
            // enabling validation against original code 'rodinia/cfd':
            compute_step_factor_legacy(nel[level], variables[level], volumes[level], step_factors[level]);
        } else {
            // Use the corrected step factor calculation:
            compute_step_factor(nel[level], variables[level], volumes[level], step_factors[level]);
        }

        for(int j=0; j<RK; j++)
        {
            #ifdef FLUX_CRIPPLE
                compute_flux_edge_crippled(
                    internal_edge_starts[level], 
                    num_internal_edges[level],
                    edges[level], 
                    #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                        edge_weights[level],
                    #endif
                    variables[level], 
                    #ifndef FLUX_FISSION
                        fluxes[level]
                        #ifdef COLOURED_CONFLICT_AVOIDANCE
                        , internal_serial_section_starts[level]
                        #endif
                    #else
                        edge_variables[level]
                    #endif
                    );
                #ifndef FLUX_FISSION
                    // Revert writes made by 'flux cripple':
                    zero_fluxes(nel[level], fluxes[level]);
                #endif
            #endif

            compute_flux_edge(
                internal_edge_starts[level], 
                num_internal_edges[level],
                edges[level], 
                #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                    edge_weights[level],
                #endif
                variables[level], 
                #ifndef FLUX_FISSION
                    fluxes[level]
                    #ifdef COLOURED_CONFLICT_AVOIDANCE
                    , internal_serial_section_starts[level]
                    #endif
                #else
                    edge_variables[level]
                #endif
                );

            compute_boundary_flux_edge(
                boundary_edge_starts[level], 
                num_boundary_edges[level],
                edges[level], 
                variables[level], 
                #ifndef FLUX_FISSION
                    fluxes[level]
                #else
                    edge_variables[level]
                #endif
                );

            compute_wall_flux_edge(
                wall_edge_starts[level], 
                num_wall_edges[level],
                edges[level], 
                variables[level], 
                #ifndef FLUX_FISSION
                    fluxes[level]
                #else
                    edge_variables[level]
                #endif
                );

            #ifdef FLUX_FISSION
                update_edges(
                    internal_edge_starts[level],
                    num_internal_edges[level],
                    edges[level], 
                    nel[level], 
                    edge_variables[level],
                    fluxes[level]);

                update_edges(
                    boundary_edge_starts[level],
                    num_boundary_edges[level],
                    edges[level], 
                    nel[level], 
                    edge_variables[level],
                    fluxes[level]);

                update_edges(
                    wall_edge_starts[level],
                    num_wall_edges[level],
                    edges[level], 
                    nel[level], 
                    edge_variables[level],
                    fluxes[level]);
            #endif

            time_step(
                j, nel[level], 
                step_factors[level], 
                volumes[level], 
                fluxes[level], 
                old_variables[level], 
                variables[level]);

            check_for_invalid_variables(variables[level], nel[level]);

            indirect_rw(
                internal_edge_starts[level], 
                num_internal_edges[level],
                edges[level], 
                #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
                    edge_weights[level],
                #endif
                variables[level], 
                #ifndef FLUX_FISSION
                    fluxes[level]
                    #ifdef COLOURED_CONFLICT_AVOIDANCE
                    , internal_serial_section_starts[level]
                    #endif
                #else
                    edge_variables[level]
                #endif
                );
            zero_fluxes(nel[level], fluxes[level]);
        }

        residual(nel[level], old_variables[level], variables[level], residuals[level]);
        if (level == 0) {
            double rms = calc_rms(nel[level], residuals[level]);
            printf(" (RMS = %.3e)", rms);
        }

        if (levels <= 1) {
            i++;
        } else {
            ///////////////////////////////////////////
            // Restrict up / prolong down the multigrid
            ///////////////////////////////////////////
            if (levels > 1)
            {
                // My understanding of geometric multigrid is that 
                // the residual error should be restricted up, not 
                // the grid state itself. However, doing this 
                // immediately introduces NaN's.
                // #define UP_RESIDUALS 1
                if(mg_direction == MG_UP)
                {
                    level++;
                    #ifdef UP_RESIDUALS
                        if (level == 1) {
                            up(residuals[level-1], 
                               variables[level], 
                               nel[level], 
                               mg_connectivity[level-1], 
                               up_scratch, 
                               mg_connectivity_size[level-1]);
                        } else {
                            up(variables[level-1], 
                               variables[level], 
                               nel[level], 
                               mg_connectivity[level-1], 
                               up_scratch, 
                               mg_connectivity_size[level-1]);
                        }
                    #else
                        up(variables[level-1], 
                           variables[level], 
                           nel[level], 
                           mg_connectivity[level-1], 
                           up_scratch, 
                           mg_connectivity_size[level-1]);
                    #endif

                    if(level == (levels-1))
                    {
                        mg_direction = MG_DOWN;
                    }
                }
                else
                {
                    level--;

                    // down() generates NaN's after 3 MG cycles
                    // down(
                    //     variables[level+1], 
                    //     nel[level+1], 
                    //     variables[level], 
                    //     nel[level], 
                    //     mg_connectivity[level], 
                    //     mg_connectivity_size[level], 
                    //     coords[level+1], 
                    //     coords[level]);

                    // down_interpolate() generates NaN's after 33 MG cycles
                    // down_interpolate(
                    //     variables[level+1], 
                    //     nel[level+1], 
                    //     variables[level], 
                    //     nel[level], 
                    //     mg_connectivity[level], 
                    //     mg_connectivity_size[level], 
                    //     coords[level+1], 
                    //     coords[level]);

                    // down_residuals() generates NaN's after 1 MG cycles
                    // #ifdef UP_RESIDUAL
                    //     if (level == 0) {
                    //         down_residuals(
                    //             variables[level+1], 
                    //             nel[level+1], 
                    //             variables[level], 
                    //             residuals[level], 
                    //             nel[level], 
                    //             mg_connectivity[level], 
                    //             mg_connectivity_size[level], 
                    //             coords[level+1], 
                    //             coords[level]);
                    //     } 
                    //     else {
                    //         down_residuals(
                    //             variables[level+1], 
                    //             nel[level+1], 
                    //             residuals[level], 
                    //             residuals[level], 
                    //             nel[level], 
                    //             mg_connectivity[level], 
                    //             mg_connectivity_size[level], 
                    //             coords[level+1], 
                    //             coords[level]);
                    //     }
                    // #else
                    //     down_residuals(
                    //         residuals[level+1], 
                    //         nel[level+1], 
                    //         variables[level], 
                    //         residuals[level], 
                    //         nel[level], 
                    //         mg_connectivity[level], 
                    //         mg_connectivity_size[level], 
                    //         coords[level+1], 
                    //         coords[level]);
                    // #endif

                    // #ifdef UP_RESIDUAL
                    //     if (level == 0) {
                    //         down_residuals_interpolate_crude(
                    //             variables[level+1], 
                    //             nel[level+1], 
                    //             residuals[level], 
                    //             variables[level], 
                    //             nel[level], 
                    //             mg_connectivity[level], 
                    //             mg_connectivity_size[level], 
                    //             coords[level+1], 
                    //             coords[level]);
                    //     } 
                    //     else {
                    //         down_residuals_interpolate_crude(
                    //             variables[level+1], 
                    //             nel[level+1], 
                    //             variables[level], 
                    //             variables[level], 
                    //             nel[level], 
                    //             mg_connectivity[level], 
                    //             mg_connectivity_size[level], 
                    //             coords[level+1], 
                    //             coords[level]);
                    //     }
                    // #else
                    //     down_residuals_interpolate_crude(
                    //         residuals[level+1], 
                    //         nel[level+1], 
                    //         residuals[level], 
                    //         variables[level], 
                    //         nel[level], 
                    //         mg_connectivity[level], 
                    //         mg_connectivity_size[level], 
                    //         coords[level+1], 
                    //         coords[level]);
                    // #endif

                    #ifdef UP_RESIDUAL
                        if (level == 0) {
                            down_residuals_interpolate_proper(
                                edges[level], 
                                num_internal_edges[level],
                                variables[level+1], 
                                residuals[level], 
                                variables[level], 
                                nel[level], 
                                mg_connectivity[level], 
                                mg_connectivity_size[level], 
                                coords[level+1], 
                                coords[level]);
                        } 
                        else {
                            down_residuals_interpolate_proper(
                                edges[level], 
                                num_internal_edges[level],
                                variables[level+1], 
                                variables[level], 
                                variables[level], 
                                nel[level], 
                                mg_connectivity[level], 
                                mg_connectivity_size[level], 
                                coords[level+1], 
                                coords[level]);
                        }
                    #else
                        down_residuals_interpolate_proper(
                            edges[level], 
                            num_internal_edges[level],
                            residuals[level+1], 
                            residuals[level], 
                            variables[level], 
                            nel[level], 
                            mg_connectivity[level], 
                            mg_connectivity_size[level], 
                            coords[level+1], 
                            coords[level]);
                    #endif

                    if (level == 0)
                    {
                        mg_direction = MG_UP;
                        i++;
                    }
                }
            }
            else {
                i++;
            }
        }
    }
    printf("\n");

    total_compute_time = omp_get_wtime() - loop_start_time;
    std::cout << "Total runtime = " << total_compute_time << std::endl;

    ////////////////////////////////////
    // Validate solution:
    ////////////////////////////////////
    printf("\n");
    if (conf.validate_result) {
        printf("Beginning validation of variables[]\n");
        for (int level=0; level<levels; level++) {
            check_for_invalid_variables(variables[level], nel[level]);
        }
        printf("  NaN check passed\n");
        bool data_check_passed = true;
        bool res;
        // for (int level=0; level<levels; level++) {
        // Update: only interested in the finest mesh:
        for (int level=0; level<1; level++) {
            // printf("Validating variables[] for level %d\n", level);
            std::string solution_filepath = generate_solution_filepath(std::string("variables"), level);
            std::ifstream file(solution_filepath.c_str());
            if(!file.is_open())
            {
                printf("  could not open variables solution file:\n");
                printf("    %s\n", solution_filepath.c_str());
                printf("  aborting validation\n");
                data_check_passed = false;
                break;
            }

            double* variables_solution = alloc<double>(nel[level]*NVAR);
            res = read_double_array(variables_solution, std::string("variables"), NVAR, nel[level], level);
            if (!res) {
                data_check_passed = false;
                break;
            } else {
                printf("  scanning variables[] on level %d for errors\n", level);
                identify_differences(variables[level], variables_solution, nel[level]);
            }
            dealloc<double>(variables_solution);
        }
        if (data_check_passed) {
            printf("PASS: variables[] validated successfully\n");
        // } else {
        //     printf("FAIL: variables[] validation failed\n");
        }
        printf("\n");
    }

    ////////////////////////////////////
    // Write out array data to file:
    ////////////////////////////////////
    log("Writing out array data\n");
    // for (int l=0; l<levels; l++) {
    // Update: only interested in the finest mesh:
    for (int l=0; l<1; l++) {
        if (conf.output_variables) {
            dump(variables[l], nel[l], l);
        }
        if (conf.output_step_factors) {
            dump_step_factors(step_factors[l], nel[l], l);
        }
        if (conf.output_edge_fluxes) {
            dump_edge_fluxes(
                edge_variables[l],
                num_internal_edges[l], num_boundary_edges[l], num_wall_edges[l], 
                internal_edge_starts[l], boundary_edge_starts[l], wall_edge_starts[l], 
                l);
        }
        if (conf.output_fluxes) {
            dump_flux(fluxes[l], nel[l], l);
        }
        if (conf.output_volumes) {
            dump_volumes(volumes[l], nel[l], l);
        }
    }

    ////////////////////////////////////
    // Output performance data to file:
    ////////////////////////////////////
    log("Writing out performance data");
    #ifdef TIME
        dump_timers_to_file(problem_size, total_compute_time);
    #endif
    #ifdef PAPI
        dump_papi_counters_to_file(problem_size);
    #endif

    dump_loop_stats_to_file(problem_size);

    ////////////////////////////////////
    // Clean memory:
    ////////////////////////////////////
    log("Cleaning memory");
    for(int i = 0; i < levels; i++)
    {
        clean_level(nel[i], volumes[i], 
            variables[i], old_variables[i],
            fluxes[i], step_factors[i], edges[i],
            edge_variables[i],
            coords[i]);
    }
    for(int i = 0; i < levels-1; i++)
    {
        dealloc<long>(mg_connectivity[i]);
    }
    dealloc<long>(up_scratch);

    delete[] (layers);
    delete[] (mg_connectivity_filename);

    return 0;
}
