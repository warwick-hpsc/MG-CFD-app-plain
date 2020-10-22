#include <errno.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef LEGACY_ORDERING
#include <algorithm>
#endif

#ifdef PAPI
#include <papi.h>
#include <set>
#endif

#include "io_enhanced.h"

std::string generate_input_binary_filename_suffix() {
    char hostname[256];
    gethostname(hostname, 256);
    std::string suffix = std::string(hostname) + std::string(".bin");
    return suffix;
}

std::string generate_output_filename_suffix(int level) {
    std::string suffix = "size=" + number_to_string(conf.mesh_duplicate_count) + "x";
    suffix += ".cycles=" + number_to_string(conf.num_cycles);
    if (level >= 0) {
        suffix += ".level=" + number_to_string(level);
    }

    return suffix;
}

std::string generate_output_filepath(std::string filename, int level)
{
    std::string filepath;
    std::string prefix(conf.output_file_prefix);
    if (prefix != "") {
        filepath += prefix;
        if (prefix.at(prefix.length()-1) != '/') {
            filepath += ".";
        }
    }
    filepath += filename;
    std::string suffix = generate_output_filename_suffix(level);
    if (suffix != "") {
        filepath += "." + suffix;
    }
    return filepath;
}

std::string generate_solution_filepath(std::string filename, int level)
{
    std::string filepath;

    std::string prefix(conf.input_file_directory);
    if (prefix.length() > 0 && prefix.at(prefix.length()-1) != '/') {
        prefix += "/";
    }
    filepath += prefix;

    filepath += "solution.";

    filepath += filename;

    std::string suffix = generate_output_filename_suffix(level);
    if (suffix != "") {
        filepath += "." + suffix;
    }

    return filepath;
}

void copy_and_shift_edges(
    edge_neighbour* source, 
    edge_neighbour* dest, 
    long num_edges,
    long node_idx_shift)
{
    for (long i=0; i<num_edges; i++) {
        dest[i] = source[i];
        if (dest[i].a >= 0) dest[i].a += node_idx_shift;
        if (dest[i].b >= 0) dest[i].b += node_idx_shift;
    }
}

void duplicate_mesh(
    long* nel,
    double** volumes,
    double3** coords,
    long* num_edges,
    long* num_iedges, 
    long* num_bedges, 
    long* num_wedges, 
    long* bedges_start, 
    long* wedges_start,
    #if defined BIN_COLOURED_VECTORS || defined BIN_COLOURED_CONTIGUOUS
        long* internal_serial_section_start,
        long* boundary_serial_section_start, 
        long* wall_serial_section_start,
    #endif
    edge_neighbour** edges,
    long nel_above,
    long** mg_mapping,
    long* mgc)
{
    log("Entering duplicate_mesh()");

    const int m = conf.mesh_duplicate_count;

    double* volumes_duplicated = alloc<double>((*nel)*m);
    double3* coords_duplicated = alloc<double3>((*nel)*m);
    #if defined OMP
    #pragma omp parallel for
    #endif
    for (int i=0; i<m; i++) {
        for (long p=0; p<(*nel); p++) {
            volumes_duplicated[ i*(*nel) + p ] = (*volumes)[p];
        }
        for (long p=0; p<(*nel); p++) {
            coords_duplicated[ i*(*nel) + p ] = (*coords)[p];
        }
    }

    long num_edges_duplicated = m*(*num_edges);

    long num_iedges_duplicated = m*(*num_iedges);
    long num_bedges_duplicated = m*(*num_bedges);
    long num_wedges_duplicated     = m*(*num_wedges);

    long bedges_start_duplicated = (*bedges_start)*m;
    long wedges_start_duplicated     = (*wedges_start)*m;

    #if defined BIN_COLOURED_VECTORS || defined BIN_COLOURED_CONTIGUOUS
        int internal_serial_section_start_duplicated = m*(*internal_serial_section_start);

        int boundary_serial_section_start_duplicated = bedges_start_duplicated + 
                m*((*boundary_serial_section_start)-(*bedges_start));

        int wall_serial_section_start_duplicated = wedges_start_duplicated + 
                m*((*wall_serial_section_start)-(*wedges_start));
    #endif

    edge_neighbour* edges_duplicated = alloc<edge_neighbour>(num_edges_duplicated);
    long j=0;
    #if defined BIN_COLOURED_VECTORS || defined BIN_COLOURED_CONTIGUOUS
        for (int i=0; i<m; i++) {
            const long n = *internal_serial_section_start;
            copy_and_shift_edges(*edges,
                                 edges_duplicated+j, 
                                 n, 
                                 (*nel)*i);
            j += n;
        }
        for (int i=0; i<m; i++) {
            const long n = (*num_iedges)-(*internal_serial_section_start);
            copy_and_shift_edges(*edges+(*internal_serial_section_start), 
                                 edges_duplicated+j, 
                                 n, 
                                 (*nel)*i);
            j += n;
        }
    #else
        for (int i=0; i<m; i++) {
            const long n = *num_iedges;
            copy_and_shift_edges(*edges,
                                 edges_duplicated+j, 
                                 n, 
                                 (*nel)*i);
            j += n;
        }
    #endif
    for (; j<bedges_start_duplicated; j++) {
        edges_duplicated[j].a = -5; 
        edges_duplicated[j].b = -5;
    }
    #if defined BIN_COLOURED_VECTORS || defined BIN_COLOURED_CONTIGUOUS
        for (int i=0; i<m; i++) {
            const long n = *boundary_serial_section_start-(*bedges_start);
            copy_and_shift_edges(*edges+(*bedges_start), 
                                 edges_duplicated+j, 
                                 n, 
                                 (*nel)*i);
            j += n;
        }
        for (int i=0; i<m; i++) {
            const long n = (*bedges_start) + (*num_bedges) - (*boundary_serial_section_start);
            copy_and_shift_edges(*edges+(*boundary_serial_section_start), 
                                 edges_duplicated+j, 
                                 n, 
                                 (*nel)*i);
            j += n;
        }
    #else
        for (int i=0; i<m; i++) {
            const long n = *num_bedges;
            copy_and_shift_edges(*edges+(*bedges_start), 
                                 edges_duplicated+j, 
                                 n, 
                                 (*nel)*i);
            j += n;
        }
    #endif
    for (; j<wedges_start_duplicated; j++) {
        edges_duplicated[j].a = -5;
        edges_duplicated[j].b = -5;
    }
    #if defined BIN_COLOURED_VECTORS || defined BIN_COLOURED_CONTIGUOUS
        for (int i=0; i<m; i++) {
            const long n = *wall_serial_section_start-(*wedges_start);
            copy_and_shift_edges(*edges+(*wedges_start), 
                                 edges_duplicated+j, 
                                 n, 
                                 (*nel)*i);
            j += n;
        }
        for (int i=0; i<m; i++) {
            const long n = (*wedges_start) + (*num_wedges) - (*wall_serial_section_start);
            copy_and_shift_edges(*edges+(*wall_serial_section_start), 
                                 edges_duplicated+j, 
                                 n, 
                                 (*nel)*i);
            j += n;
        }
    #else
        for (int i=0; i<m; i++) {
            const long n = *num_wedges;
            copy_and_shift_edges(*edges+(*wedges_start), 
                                 edges_duplicated+j, 
                                 n, 
                                 (*nel)*i);
            j += n;
        }
    #endif
    for (; j<num_edges_duplicated; j++) {
        edges_duplicated[j].a = -5;
        edges_duplicated[j].b = -5;
    }

    if (mg_mapping != NULL) {
        long mgc_duplicated = m*(*mgc);
        long* mg_mapping_duplicated = alloc<long>(mgc_duplicated);
        for (int i=0; i<m; i++) {
            for (long n=0; n<(*mgc); n++) {
                mg_mapping_duplicated[i*(*mgc) + n] = (*mg_mapping)[n] + (nel_above*i);
            }
        }
        dealloc<long>(*mg_mapping);
        *mg_mapping = mg_mapping_duplicated;
        *mgc = mgc_duplicated;
    }

    *nel *= m;

    dealloc<double>(*volumes);
    *volumes = volumes_duplicated;
    dealloc<double3>(*coords);
    *coords = coords_duplicated;

    *num_edges = num_edges_duplicated;
    *num_iedges = num_iedges_duplicated;
    *num_bedges = num_bedges_duplicated;
    *num_wedges     = num_wedges_duplicated;
    *bedges_start = bedges_start_duplicated;
    *wedges_start     = wedges_start_duplicated;

    dealloc<edge_neighbour>(*edges);
    *edges = edges_duplicated;

    #if defined BIN_COLOURED_VECTORS || defined BIN_COLOURED_CONTIGUOUS
        *internal_serial_section_start = internal_serial_section_start_duplicated;
        *boundary_serial_section_start = boundary_serial_section_start_duplicated;
        *wall_serial_section_start     = wall_serial_section_start_duplicated;
    #endif

    log("Exiting inflate_mesh()");
}

bool read_grid_from_bin(
    const char* data_file_name, 
    long* nel, 
    double** volumes, 
    long* num_edges, 
    #ifdef E_SOA
        long* edges_soa_step, 
    #endif
    long* num_iedges, 
    long* num_bedges, 
    long* num_wedges, 
    long* iedges_start, 
    long* bedges_start, 
    long* wedges_start, 
    edge_neighbour** edges, 
    double3** coords, 
    long** mg_connectivity, 
    long* mg_size)
{
    log("read_grid_from_bin() called");
    log("Attempting to open '%s'", data_file_name);
    FILE* fp;
    if ((fp=fopen(data_file_name, "rb")) == NULL) {
        log("'%s' binary file cannot be read", data_file_name);
        return false;
    }

    if (fread(nel, sizeof(long), 1, fp) != 1) {
        log("Corruption detected in '%s'", data_file_name);
        fclose(fp);
        return false;
    }
    if (fread(num_edges, sizeof(long), 1, fp) != 1) {
        log("Corruption detected in '%s'", data_file_name);
        fclose(fp);
        return false;
    }
    if (fread(num_iedges, sizeof(long), 1, fp) != 1) {
        log("Corruption detected in '%s'", data_file_name);
        fclose(fp);
        return false;
    }
    if (fread(num_bedges, sizeof(long), 1, fp) != 1) {
        log("Corruption detected in '%s'", data_file_name);
        fclose(fp);
        return false;
    }
    if (fread(num_wedges, sizeof(long), 1, fp) != 1) {
        log("Corruption detected in '%s'", data_file_name);
        fclose(fp);
        return false;
    }
    if (fread(iedges_start, sizeof(long), 1, fp) != 1) {
        log("Corruption detected in '%s'", data_file_name);
        fclose(fp);
        return false;
    }
    if (fread(bedges_start, sizeof(long), 1, fp) != 1) {
        log("Corruption detected in '%s'", data_file_name);
        fclose(fp);
        return false;
    }
    if (fread(wedges_start, sizeof(long), 1, fp) != 1) {
        log("Corruption detected in '%s'", data_file_name);
        fclose(fp);
        return false;
    }

    if ((*nel) < 0) {
        log("Corruption detected in '%s'", data_file_name);
        return false;
    }
    unsigned long nel_u = (unsigned long)(*nel);

    // Check that the above variables were loaded correctly. Any 
    // inconsistency is indicative of the bin file having been 
    // created on a different architecture that encodes differently 
    // and so cannot be loaded here:
    if (*nel < 0 || *num_edges < 0 || *num_iedges < 0 ||
        *num_bedges < 0 || *num_wedges < 0 ||
        *iedges_start < 0 || *bedges_start < 0 || 
        *wedges_start < 0 || 
        *num_iedges + *num_bedges + *num_wedges > *num_edges || 
        *iedges_start > *num_edges || 
        *bedges_start > *num_edges || 
        *wedges_start > *num_edges) {
        log("Corruption detected in '%s'", data_file_name);
        fclose(fp);
        return false;
    }
    unsigned long num_edges_u = (unsigned long)(*num_edges);

    *volumes = alloc<double>(nel_u);
    if (fread(*volumes, sizeof(double), *nel, fp) != nel_u) {
        log("Corruption detected in '%s'", data_file_name);
        fclose(fp);
        return false;
    }

    *edges = alloc<edge_neighbour>(*num_edges);
    if (fread(*edges, sizeof(edge_neighbour), *num_edges, fp) != num_edges_u) {
        log("Corruption detected in '%s'", data_file_name);
        fclose(fp);
        return false;
    }

    #ifdef LEGACY_ORDERING
        std::sort(*edges, 
                  (*edges)+(*num_iedges), 
                  compare_two_edges);
        std::sort(*edges+(*num_iedges), 
                  (*edges)+(*num_iedges)+(*num_bedges), 
                  compare_two_edges);
        std::sort(*edges+(*num_iedges)+(*num_bedges), 
                  (*edges)+(*num_iedges)+(*num_bedges)+(*num_wedges), 
                  compare_two_edges);
    #endif

    *coords = alloc<double3>(nel_u);
    if (fread(*coords, sizeof(double3), *nel , fp) != nel_u) {
        log("Corruption detected in '%s'", data_file_name);
        fclose(fp);
        return false;
    }

    if (fread(mg_size, sizeof(long), 1, fp) != 1) {
        log("Corruption detected in '%s'", data_file_name);
        fclose(fp);
        return false;
    }
    if ((*mg_size) < 0) {
        *mg_size = 0;
        *mg_connectivity = NULL;
        log("'%s' cannot be read\n", data_file_name);
        fclose(fp);
        return false;
    }
    else if (*mg_size != 0) {
        unsigned long mg_size_u = (unsigned long)(*mg_size);
        *mg_connectivity = alloc<long>(mg_size_u);
        if (fread(*mg_connectivity, sizeof(long), *mg_size, fp) != mg_size_u) {
            log("Corruption detected in '%s'", data_file_name);
            fclose(fp);
            return false;
        }
    } else {
        *mg_connectivity = NULL;
    }

    fclose(fp);

    return true;
}

bool write_grid_to_bin(
    const char* data_file_name, 
    long nel, 
    const double *restrict volumes, 
    long num_edges, 
    long num_iedges, 
    long num_bedges, 
    long num_wedges, 
    long iedges_start, 
    long bedges_start, 
    long wedges_start, 
    const edge_neighbour* edges, 
    const double3* coords, 
    const long* mg_connectivity, 
    long mg_size)
{
    log("write_grid_to_bin() called");

    log("Writing to '%s'\n", data_file_name);

    FILE* fp;
    fp = fopen(data_file_name, "wb"); 
    if (fp == NULL) {
        fprintf(stderr, "ERROR: Could not open '%s' for writing: %s\n", data_file_name, strerror(errno));
        DEBUGGABLE_ABORT
    }

    fwrite(&nel, sizeof(long), 1, fp);
    fwrite(&num_edges, sizeof(long), 1, fp);
    fwrite(&num_iedges, sizeof(long), 1, fp);
    fwrite(&num_bedges, sizeof(long), 1, fp);
    fwrite(&num_wedges, sizeof(long), 1, fp);
    fwrite(&iedges_start, sizeof(long), 1, fp);
    fwrite(&bedges_start, sizeof(long), 1, fp);
    fwrite(&wedges_start, sizeof(long), 1, fp);

    fwrite(volumes, sizeof(double), nel, fp);
    fwrite(edges, sizeof(edge_neighbour), num_edges, fp);
    fwrite(coords, sizeof(double3), nel, fp);

    fwrite(&mg_size, sizeof(long), 1, fp);
    if (mg_size > 0) {
        fwrite(mg_connectivity, sizeof(long), mg_size, fp);
    }

    fclose(fp);

    return true;
}

void read_input_dat(
    const char* file_name, 
    int* problem_size, 
    std::string** layers, 
    std::string** mg_connectivity)
{
    std::ifstream file(file_name);

    bool have_size=false, have_num_levels=false, have_mesh_name=false, have_mesh_filenames=false, have_mg_mapping=false;
    std::string file_line;

    if (!file.is_open()) {
        fprintf(stderr, "Error: Could not open input file '%s'\n", file_name);
        DEBUGGABLE_ABORT
    }
    else {
        while (std::getline(file, file_line))
        {
            if (file_line.c_str()[0] == '#')
                continue;

            std::istringstream str_iss(file_line);
            if (file_line.c_str()[0] == '[')
            {
                if (strcmp(file_line.c_str(), "[levels]")==0)
                {
                    if (!have_num_levels)
                    {
                        fprintf(stderr, "Error parsing %s: Need to know number of levels before parsing level filenames\n", file_name);
                        DEBUGGABLE_ABORT
                    }
                    have_mesh_filenames = true;
                    *layers = new std::string[levels];
                    for (int i=0; i<levels; i++)
                    {
                        std::string key;
                        std::string value;
                        if (!std::getline(file, file_line))
                        {
                            fprintf(stderr, "Error parsing %s: Have reached EOF before reading all mesh filenames\n", file_name);
                            DEBUGGABLE_ABORT
                        }
                        else
                        {
                            str_iss.clear();
                            str_iss.seekg(0, std::ios::beg);
                            str_iss.str(file_line);
                            if (!std::getline(str_iss, key, '='))
                            {
                                fprintf(stderr, "Error parsing '%s': Was expecting a key-value pair following [levels]\n", file_name);
                                DEBUGGABLE_ABORT
                            }
                            else
                            {
                                if (std::getline(str_iss, value)) {
                                    key = trim(key);
                                    value = trim(value);
                                    int idx = atoi(key.c_str());
                                    (*layers)[idx] = strdup(value.c_str());
                                }
                            }
                        }
                    }
                }
                else if (strcmp(file_line.c_str(), "[mg_mapping]")==0)
                {
                    if (!have_num_levels) {
                        fprintf(stderr, "Error parsing %s: Need to know number of levels before parsing level filenames\n", file_name);
                        DEBUGGABLE_ABORT
                    }
                    have_mg_mapping = true;
                    *mg_connectivity = new std::string[levels-1];
                    for (int i=0; i<levels-1; i++)
                    {
                        std::string key;
                        std::string value;
                        if (!std::getline(file, file_line)) {
                            fprintf(stderr, "Error parsing %s: Have reached EOF before reading all MG filenames\n", file_name);
                            DEBUGGABLE_ABORT
                        }
                        else
                        {
                            str_iss.clear();
                            str_iss.seekg(0, std::ios::beg);
                            str_iss.str(file_line);
                            if (!std::getline(str_iss, key, '=')) {
                                fprintf(stderr, "Error parsing '%s': Was expecting a key-value pair following [mg_mapping]\n", file_name);
                                DEBUGGABLE_ABORT
                            } else {
                                if (std::getline(str_iss, value)) {
                                    key = trim(key);
                                    value = trim(value);
                                    int idx = atoi(key.c_str());
                                    (*mg_connectivity)[idx] = strdup(value.c_str());
                                }
                            }
                        }
                    }
                    have_mg_mapping = true;
                }
            }
            else
            {
                std::string key;
                std::string value;
                if (std::getline(str_iss, key, '=')) {
                    // Line contains a key-value pair
                    if (std::getline(str_iss, value))
                    {
                        key = trim(key);
                        value = trim(value);

                        if (strcmp(key.c_str(), "size")==0) {
                            *problem_size = atoi(value.c_str());
                            have_size = true;
                        }
                        else if (strcmp(key.c_str(), "num_levels")==0) {
                            levels = atoi(value.c_str());
                            have_num_levels = true;
                        }
                        else if (strcmp(key.c_str(), "mesh_name")==0) {
                            if (strcmp(value.c_str(), "la_cascade")==0) {
                                mesh_variant = MESH_LA_CASCADE;
                                have_mesh_name = true;
                            }
                            else if (strcmp(value.c_str(), "rotor37")==0) {
                                mesh_variant = MESH_ROTOR_37;
                                have_mesh_name = true;
                            }
                            else if (strcmp(value.c_str(), "fvcorr")==0) {
                                mesh_variant = MESH_FVCORR;
                                have_mesh_name = true;
                            }
                            else if (strcmp(value.c_str(), "m6wing")==0) {
                                mesh_variant = MESH_M6_WING;
                                have_mesh_name = true;
                            }
                            else {
                                fprintf(stderr, "Error parsing %s: Unknown mesh_name '%s'\n", file_name, value.c_str());
                                DEBUGGABLE_ABORT
                            }
                        }
                    }
                }
            }
        }
    }

    if (!have_size) {
        fprintf(stderr, "Error parsing '%s': size not present\n", file_name);
        DEBUGGABLE_ABORT
    }
    if (!have_num_levels) {
        fprintf(stderr, "Error parsing '%s': number of levels not present\n", file_name);
        DEBUGGABLE_ABORT
    }
    if (!have_mesh_name) {
        fprintf(stderr, "Error parsing '%s': mesh name not present\n", file_name);
        DEBUGGABLE_ABORT
    }
    if (!have_mesh_filenames) {
        fprintf(stderr, "Error parsing '%s': mesh filenames not present\n", file_name);
        DEBUGGABLE_ABORT
    }

    if (!have_mg_mapping) {
        *mg_connectivity = new std::string[levels-1];
        for (int i=0; i<levels-1; i++) {
            (*mg_connectivity)[i] = (char*)malloc(sizeof(char));
            (*mg_connectivity)[i][0] = '\0';
        }
    }
}

#ifdef PAPI
void read_papi_config(const char* filename, int* nevents, int** events)
{
    if (access(filename, F_OK) == -1) {
        fprintf(stderr, "ERROR: PAPI filepath \"%s\" does not exist.\n", filename);
        fprintf(stderr, "%s:%d\n", __FILE__, __LINE__);
        DEBUGGABLE_ABORT
    }

    int code, ret;
    std::set<int> events_set;
    std::ifstream file_if(filename);
    std::string line;
    std::string events_concat = "";
    *nevents = 0;
    while(std::getline(file_if, line))
    {
        if (line.c_str()[0] == '#' || strcmp(line.c_str(), "")==0) {
            continue;
        }

        ret = PAPI_event_name_to_code(strdup(line.c_str()), &code);
        if (ret != PAPI_OK) {
            printf("Could not convert string '%s' to PAPI event, error = %s\n", line.c_str(), PAPI_strerror(ret));
        } else {
            if (PAPI_query_event(code) != PAPI_OK) {
                printf("PAPI event %s not present\n", line.c_str());
            } else {
                events_set.insert(code);
                if (strcmp(events_concat.c_str(), "")==0) {
                    events_concat += line;
                } else {
                    events_concat += std::string(", ") + line;
                }
            }
        }
    }
    printf("Monitoring PAPI events: %s\n", events_concat.c_str());

    *nevents = events_set.size();
    *events = alloc<int>(events_set.size());
    int i=0;
    for (std::set<int>::iterator iter=events_set.begin(); iter!=events_set.end(); iter++, i++) {
        (*events)[i] = *iter;
    }
}
#endif

void read_mg_connectivity(const char* file_name, long** mg_connectivity, long* mgc)
{
    log("Opening '%s'\n", file_name);

    std::ifstream file(file_name);
    *mgc = 0;
    if(file.is_open())
    {
        file >> *mgc;
        *mg_connectivity = alloc<long>(*mgc);

        for(long i = 0; i < *mgc; i++)
        {
            file >> (*mg_connectivity)[i];
        }
    }
    else
    {
        std::cout << "could not open mg file: '" << file_name << "'" << std::endl;
        DEBUGGABLE_ABORT
    }
}

void dump_step_factors(
    const double *restrict step_factors, 
    long npoints, int level)
{
    log("dump_step_factors() called");

    std::string filepath = generate_output_filepath(std::string("step_factors"), level);
    FILE *file = fopen(filepath.c_str(), "w");
    if (file == NULL) {
        fprintf(stderr, "ERROR: Failed to open file for writing: '%s'\n", filepath.c_str());
        exit(EXIT_FAILURE);
    }
    for(long i=0; i<npoints; i++) {
        fprintf(file, "%.17e\n", step_factors[i]);
    }
    fclose(file);

    log("dump_step_factors() complete");
}

void dump_edge_fluxes(
    const edge* edge_variables,
    long num_iedges, 
    long num_bedges, 
    long num_wedges, 
    long iedges_start, 
    long bedges_start, 
    long wedges_start, 
    int level)
{
    log("dump_edge_fluxes() called");

    if (edge_variables == NULL) {
        #ifndef FLUX_FISSION
            return;
        #else
            printf("ERROR: edge_variables is NULL\n");
            DEBUGGABLE_ABORT
        #endif
    }

    {
        std::string mx_filepath = generate_output_filepath(std::string("edge_mx"), level);
        FILE* file = fopen(mx_filepath.c_str(), "w");
        if (file == NULL) {
            fprintf(stderr, "ERROR: Failed to open file for writing: '%s'\n", mx_filepath.c_str());
            exit(EXIT_FAILURE);
        }
        for(long i=iedges_start; i<iedges_start+num_iedges; i++) {
            fprintf(file, "%.17e %.17e\n", edge_variables[i*NVAR + VAR_MOMENTUMX].a, edge_variables[i*NVAR + VAR_MOMENTUMX].b);
        }
        for(long i=bedges_start; i<bedges_start+num_bedges; i++) {
            fprintf(file, "%.17e %.17e\n", edge_variables[i*NVAR + VAR_MOMENTUMX].a, edge_variables[i*NVAR + VAR_MOMENTUMX].b);
        }
        for(long i=wedges_start; i<wedges_start+num_wedges; i++) {
            fprintf(file, "%.17e %.17e\n", edge_variables[i*NVAR + VAR_MOMENTUMX].a, edge_variables[i*NVAR + VAR_MOMENTUMX].b);
        }
        fclose(file);
    }

    {
        std::string my_filepath = generate_output_filepath(std::string("edge_my"), level);
        FILE* file = fopen(my_filepath.c_str(), "w");
        if (file == NULL) {
            fprintf(stderr, "ERROR: Failed to open file for writing: '%s'\n", my_filepath.c_str());
            exit(EXIT_FAILURE);
        }
        for(long i=iedges_start; i<iedges_start+num_iedges; i++) {
            fprintf(file, "%.17e %.17e\n", edge_variables[i*NVAR + VAR_MOMENTUMY].a, edge_variables[i*NVAR + VAR_MOMENTUMY].b);
        }
        for(long i=bedges_start; i<bedges_start+num_bedges; i++) {
            fprintf(file, "%.17e %.17e\n", edge_variables[i*NVAR + VAR_MOMENTUMY].a, edge_variables[i*NVAR + VAR_MOMENTUMY].b);
        }
        for(long i=wedges_start; i<wedges_start+num_wedges; i++) {
            fprintf(file, "%.17e %.17e\n", edge_variables[i*NVAR + VAR_MOMENTUMY].a, edge_variables[i*NVAR + VAR_MOMENTUMY].b);
        }
        fclose(file);
    }

    {
        std::string mz_filepath = generate_output_filepath(std::string("edge_mz"), level);
        FILE* file = fopen(mz_filepath.c_str(), "w");
        if (file == NULL) {
            fprintf(stderr, "ERROR: Failed to open file for writing: '%s'\n", mz_filepath.c_str());
            exit(EXIT_FAILURE);
        }
        for(long i=iedges_start; i<iedges_start+num_iedges; i++) {
            fprintf(file, "%.17e %.17e\n", edge_variables[i*NVAR + VAR_MOMENTUMZ].a, edge_variables[i*NVAR + VAR_MOMENTUMZ].b);
        }
        for(long i=bedges_start; i<bedges_start+num_bedges; i++) {
            fprintf(file, "%.17e %.17e\n", edge_variables[i*NVAR + VAR_MOMENTUMZ].a, edge_variables[i*NVAR + VAR_MOMENTUMZ].b);
        }
        for(long i=wedges_start; i<wedges_start+num_wedges; i++) {
            fprintf(file, "%.17e %.17e\n", edge_variables[i*NVAR + VAR_MOMENTUMZ].a, edge_variables[i*NVAR + VAR_MOMENTUMZ].b);
        }
        fclose(file);
    }

    {
        std::string density_filepath = generate_output_filepath(std::string("edge_p"), level);
        FILE* file = fopen(density_filepath.c_str(), "w");
        if (file == NULL) {
            fprintf(stderr, "ERROR: Failed to open file for writing: '%s'\n", density_filepath.c_str());
            exit(EXIT_FAILURE);
        }
        for(long i=iedges_start; i<iedges_start+num_iedges; i++) {
            fprintf(file, "%.17e %.17e\n", edge_variables[i*NVAR + VAR_DENSITY].a, edge_variables[i*NVAR + VAR_DENSITY].b);
        }
        for(long i=bedges_start; i<bedges_start+num_bedges; i++) {
            fprintf(file, "%.17e %.17e\n", edge_variables[i*NVAR + VAR_DENSITY].a, edge_variables[i*NVAR + VAR_DENSITY].b);
        }
        for(long i=wedges_start; i<wedges_start+num_wedges; i++) {
            fprintf(file, "%.17e %.17e\n", edge_variables[i*NVAR + VAR_DENSITY].a, edge_variables[i*NVAR + VAR_DENSITY].b);
        }
        fclose(file);
    }

    {
        std::string density_energy_filepath = generate_output_filepath(std::string("edge_pe"), level);
        FILE* file = fopen(density_energy_filepath.c_str(), "w");
        if (file == NULL) {
            fprintf(stderr, "ERROR: Failed to open file for writing: '%s'\n", density_energy_filepath.c_str());
            exit(EXIT_FAILURE);
        }
        for(long i=iedges_start; i<iedges_start+num_iedges; i++) {
            fprintf(file, "%.17e %.17e\n", edge_variables[i*NVAR + VAR_DENSITY_ENERGY].a, edge_variables[i*NVAR + VAR_DENSITY_ENERGY].b);
        }
        for(long i=bedges_start; i<bedges_start+num_bedges; i++) {
            fprintf(file, "%.17e %.17e\n", edge_variables[i*NVAR + VAR_DENSITY_ENERGY].a, edge_variables[i*NVAR + VAR_DENSITY_ENERGY].b);
        }
        for(long i=wedges_start; i<wedges_start+num_wedges; i++) {
            fprintf(file, "%.17e %.17e\n", edge_variables[i*NVAR + VAR_DENSITY_ENERGY].a, edge_variables[i*NVAR + VAR_DENSITY_ENERGY].b);
        }
        fclose(file);
    }

    log("dump_edge_fluxes() complete");
}

void dump_flux(
    const double *restrict fluxes, 
    long nel, int level)
{
    log("dump_flux() called");

    std::string filepath = generate_output_filepath(std::string("fluxes"), level);
    FILE* file = fopen(filepath.c_str(), "w");
    if (file == NULL) {
        fprintf(stderr, "ERROR: Failed to open file for writing: '%s'\n", filepath.c_str());
        exit(EXIT_FAILURE);
    }

    long chunk_size = nel;
    for (long i=0; i<chunk_size; i++) {
        const long p_idx = i*NVAR + VAR_DENSITY;
        const long pe_idx = i*NVAR + VAR_DENSITY_ENERGY;
        const long mx_idx = i*NVAR + VAR_MOMENTUMX;
        const long my_idx = i*NVAR + VAR_MOMENTUMY;
        const long mz_idx = i*NVAR + VAR_MOMENTUMZ;

        fprintf(file, "%.17e %.17e %.17e %.17e %.17e\n", fluxes[p_idx], fluxes[mx_idx], fluxes[my_idx], fluxes[mz_idx], fluxes[pe_idx]);
    }
    fclose(file);

    log("dump_flux() complete");
}

void dump_volumes(
    const double *restrict volumes, 
    long npoints, int level)
{
    log("dump_volumes() called");

    std::string filepath = generate_output_filepath(std::string("volumes"), level);
    FILE *file = fopen(filepath.c_str(), "w");
    if (file == NULL) {
        fprintf(stderr, "ERROR: Failed to open file for writing: '%s'\n", filepath.c_str());
        exit(EXIT_FAILURE);
    }
    for (long i=0; i<npoints; i++) {
        fprintf(file, "%.17e\n", volumes[i]);
    }
    fclose(file);

    log("dump_volumes() complete");
}

bool read_double_array(
    double* data_array, std::string name, 
    int ndim, long nel, int level)
{
    std::string filepath = generate_solution_filepath(name, level);
    if (!file_exists(filepath.c_str())) {
        return false;
    }

    std::ifstream file(filepath.c_str());
    for (long i=0; i<nel; i++) {
        for (int v=0; v<ndim; v++) {
            file >> data_array[i*ndim + v];
        }
    }

    return true;
}

void prepare_csv_identification(
    bool write_header,
    std::string *header_s,
    std::string *data_line_s,
    int input_size)
{
    log("prepare_csv_identification() called");

    std::ostringstream header, data_line;

    if (write_header) header << "Size," ;
    data_line << input_size << "," ;

    if (write_header) header << "Mesh," ;
    if (mesh_variant == MESH_LA_CASCADE) {
        data_line << "la_cascade,";
    } else if (mesh_variant == MESH_ROTOR_37) {
        data_line << "rotor37,";
    } else if (mesh_variant == MESH_FVCORR) {
        data_line << "fvcorr,";
    } else if (mesh_variant == MESH_M6_WING) {
        data_line << "m6wing,";
    } else {
        data_line << "unknown,";
    }

    if (write_header) header << "MG cycles,";
    data_line << conf.num_cycles << "," ;

    if (write_header) header << "Flux variant,";
    #ifdef FLUX_CRIPPLE
        data_line << "FluxCripple" ;
    #else
        data_line << "Normal" ;
    #endif
    data_line << "," ;

    if (write_header) header << "Flux options,";
    #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
        data_line << "PrecomputeLength;" ;
    #endif
    #ifdef FLUX_REUSE_DIV
        data_line << "Reciprocal;" ;
    #endif
    #ifdef FLUX_REUSE_FACTOR
        data_line << "ReuseFactor;" ;
    #endif
    #ifdef FLUX_REUSE_FLUX
        data_line << "ReuseFluxes;" ;
    #endif
    data_line << "," ;

    if (write_header) header << "CC,";
    #ifdef __ICC
        data_line << "intel" << "," ;
    #elif defined __clang__
        data_line << "clang" << "," ;
    #elif defined _CRAYC
        data_line << "cray" << "," ;
    #elif defined __GNUC__
        data_line << "gnu" << "," ;
    #else
        data_line << "UNKNOWN," ;
    #endif

    if (write_header) header << "CC version,";
    #ifdef __ICC
        int intel_update_minor = __INTEL_COMPILER_UPDATE;
        if (__ICC == 1900) {
            // In version 19, __INTEL_COMPILER_UPDATE is no longer being 
            // updated to reflect the minor version. Until Intel fix this, 
            // need to figure out minor version a different way.
            if (__INTEL_COMPILER_BUILD_DATE == 20190117) {
                intel_update_minor = 2;
            } else if (__INTEL_COMPILER_BUILD_DATE == 20190206) {
                intel_update_minor = 3;
            } else if (__INTEL_COMPILER_BUILD_DATE == 20190416) {
                intel_update_minor = 4;
            } else {
                printf("WARNING: Cannot determine minor version of this Intel compiler\n");
                intel_update_minor = -1;
            }
        }
        data_line << XMACRO_TO_STR(__ICC) << "u" << intel_update_minor << "," ;
    #elif defined __clang__
        data_line << XMACRO_TO_STR(__clang_major__) << "." << XMACRO_TO_STR(__clang_minor__) << "." << XMACRO_TO_STR(__clang_patchlevel__) << "," ;
    #elif defined _CRAYC
        data_line << XMACRO_TO_STR(_RELEASE) << "." << XMACRO_TO_STR(_RELEASE_MINOR) << "," ;
    #elif defined __GNUC__
        data_line << XMACRO_TO_STR(__GNUC__) << "." << XMACRO_TO_STR(__GNUC_MINOR__) << "." << XMACRO_TO_STR(__GNUC_PATCHLEVEL__) << "," ;
    #else
        data_line << "UNKNOWN," ;
    #endif

    if (write_header) header << "Opt level," ;
    #ifdef OPT_LEVEL
        data_line << XMACRO_TO_STR(OPT_LEVEL) << "," ;
    #else
        data_line << "3" << "," ;
    #endif

    if (write_header) header << "Instruction set," ;
    #ifdef INSN_SET
        data_line << XMACRO_TO_STR(INSN_SET) << "," ;
    #else
        data_line << "UNKNOWN," ;
    #endif

    if (write_header) header << "Precise FP,";
    #ifdef PRECISE_FP
        data_line << "Y,";
    #else
        data_line << "N,";
    #endif

    if (write_header) header << "SIMD," ;
    #ifdef SIMD
        data_line << "Y,";
    #else
        data_line << "N,";
    #endif

    if (write_header) header << "SIMD conflict avoidance strategy," ;
    #ifdef COLOURED_CONFLICT_AVOIDANCE
        #ifdef BIN_COLOURED_VECTORS
            data_line << "ColouredEdgeVectors,";
        #elif defined BIN_COLOURED_CONTIGUOUS
            data_line << "ColouredEdgesContiguous,";
        #else
            data_line << "None,";
        #endif
    #elif (defined MANUAL_GATHER) && (defined MANUAL_SCATTER)
        data_line << "ManualGatherScatter,";
    #elif defined MANUAL_GATHER
        data_line << "ManualGather,";
    #elif defined MANUAL_SCATTER
        data_line << "ManualScatter,";
    #else
        data_line << "None,";
    #endif

    if (write_header) header << "SIMD len," ;
    #ifdef DBLS_PER_SIMD
        data_line << DBLS_PER_SIMD << ",";
    #else
        data_line << "1," ;
    #endif

    if (write_header) header << "OpenMP," ;
    #ifdef OMP
        data_line << "Strong," ;
    #else
        data_line << "Off," ;
    #endif

    if (write_header) header << "Num threads," ;
    #if defined OMP
        int num_threads = conf.omp_num_threads;
    #else
        int num_threads = 1;
    #endif
    data_line << number_to_string(num_threads) + ",";

    if (write_header) header << "Permit scatter OpenMP," ;
    #ifdef OMP_SCATTERS
        data_line << "Y,";
    #else
        data_line << "N,";
    #endif

    if (write_header) header << "Flux fission," ;
    #ifdef FLUX_FISSION
        data_line << "Y," ;
    #else
        data_line << "N," ;
    #endif

    if (write_header) header << "CPU," ;
    data_line << get_cpu_model_name() << "," ;

    (*header_s) = header.str();
    (*data_line_s) = data_line.str();

    log("prepare_csv_identification() complete");
}
