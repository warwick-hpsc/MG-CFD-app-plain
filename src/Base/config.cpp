//************************************************//
// Copyright 2016-2019 University of Warwick

// Permission is hereby granted, free of charge, to any person obtaining 
// a copy of this software and associated documentation files (the "Software"), 
// to deal in the Software without restriction, including without limitation 
// the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or 
// sell copies of the Software, and to permit persons to whom the Software is furnished 
// to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included 
// in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
// PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
// COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
// IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//************************************************//

#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "config.h"
#include "common.h"

config conf;

void set_config_defaults() {
    conf.config_filepath = (char*)malloc(sizeof(char));
    conf.config_filepath[0] = '\0';
    conf.input_file = (char*)malloc(sizeof(char));
    conf.input_file[0] = '\0';
    conf.input_file_directory = (char*)malloc(sizeof(char));
    conf.input_file_directory[0] = '\0';
    conf.papi_config_file = (char*)malloc(sizeof(char));
    conf.papi_config_file[0] = '\0';
    conf.output_file_prefix = (char*)malloc(sizeof(char));
    conf.output_file_prefix[0] = '\0';

    conf.mesh_duplicate_count = 1;

    conf.num_cycles = 25;

    conf.validate_result = false;

    #ifdef OMP
    conf.omp_num_threads = omp_get_max_threads();
    #else
    conf.omp_num_threads = 1;
    #endif

    conf.output_variables = false;
    conf.output_old_variables = false;
    conf.output_fluxes = false;
    conf.output_volumes = false;
    conf.output_step_factors = false;
    conf.output_edge_fluxes = false;
}

void set_config_param(const char* const key, const char* const value) {
    if (strcmp(key,"config_filepath")==0) {
        if (conf.config_filepath != NULL)
            free(conf.config_filepath);
        conf.config_filepath = strdup(value);
    }
    else if (strcmp(key,"input_file")==0) {
        if (conf.input_file != NULL)
            free(conf.input_file);
        conf.input_file = strdup(value);
    }
    else if (strcmp(key,"input_file_directory")==0) {
        if (conf.input_file_directory != NULL)
            free(conf.input_file_directory);
        conf.input_file_directory = strdup(value);
    }
    else if (strcmp(key,"papi_config_file")==0) {
        if (conf.papi_config_file != NULL)
            free(conf.papi_config_file);
        conf.papi_config_file = strdup(value);
    }
    else if (strcmp(key,"output_file_prefix")==0) {
        if (conf.output_file_prefix != NULL)
            free(conf.output_file_prefix);
        conf.output_file_prefix = strdup(value);
    }

    else if (strcmp(key, "mesh_duplicate_count")==0) {
        conf.mesh_duplicate_count = atoi(value);
    }

    else if (strcmp(key, "cycles")==0) {
        conf.num_cycles = atoi(value);
    }

    else if (strcmp(key,"omp_num_threads")==0) {
        #if defined OMP
            conf.omp_num_threads = atoi(value);
            if (conf.omp_num_threads < 1)
                conf.omp_num_threads = 1;
        #endif
    }

    else if (strcmp(key,"output_variables")==0) {
        if (strcmp(value, "Y")==0) {
            conf.output_variables = true;
        }
    }
    else if (strcmp(key,"output_old_variables")==0) {
        if (strcmp(value, "Y")==0) {
            conf.output_old_variables = true;
        }
    }
    else if (strcmp(key,"output_step_factors")==0) {
        if (strcmp(value, "Y")==0) {
            conf.output_step_factors = true;
        }
    }
    else if (strcmp(key,"output_edge_fluxes")==0) {
        if (strcmp(value, "Y")==0) {
            conf.output_edge_fluxes = true;
        }
    }
    else if (strcmp(key,"output_fluxes")==0) {
        if (strcmp(value, "Y")==0) {
            conf.output_fluxes = true;
        }
    }
    else if (strcmp(key,"output_volumes")==0) {
        if (strcmp(value, "Y")==0) {
            conf.output_volumes = true;
        }
    }
    else {
        printf("WARNING: Unknown key '%s' encountered during parsing of config file.\n", key);
    }
}

void read_config() {
    log("read_config() called\n");
    
    if (access(conf.config_filepath, F_OK) == -1) {
        fprintf(stderr, "ERROR: \"%s\" does not exist.\n", conf.config_filepath);
        // DEBUGGABLE_ABORT
        return;
    }

    std::ifstream file(conf.config_filepath);
    std::string str;
    while(std::getline(file, str))
    {
        if (str.c_str()[0] == '#')
            continue;

        log("Processing line: '%s'", str.c_str());

        std::istringstream str_iss(str);
        std::string key;
        if (std::getline(str_iss, key, '=')) {
            std::string value;
            if (std::getline(str_iss, value)) {
                key=trim(key);
                value=trim(value);
                set_config_param(key.c_str(), value.c_str());
            }
        }
    }

    #ifdef PAPI
        if (conf.papi_config_file[0] == '\0') {
            printf("ERROR: PAPI enabled but 'papi_config_file' not set in config\n");
            exit(EXIT_FAILURE);
        }
    #endif

    std::string config_dirpath;
    const size_t last_slash_idx = std::string(conf.config_filepath).rfind('/');
    // printf("last_slash_idx: %d\n", last_slash_idx);
    if (last_slash_idx != std::string::npos) {
        config_dirpath = std::string(conf.config_filepath).substr(0, last_slash_idx);
    } else {
        config_dirpath = std::string("");
    }

    if (conf.input_file_directory[0] != '/' && 
        config_dirpath != std::string("")) {
            // 'input_file_directory' is currently relative to config 
            if (std::string(conf.input_file_directory) == "./") {
                set_config_param("input_file_directory", config_dirpath.c_str());
            } else {
                // need to prepend config dirpath:
                std::string input_file_dir_corrected = std::string(config_dirpath) + "/" + conf.input_file_directory;

                set_config_param("input_file_directory", input_file_dir_corrected.c_str());
            }
    }
}

bool parse_arguments(int argc, char** argv) {
    int optc;
    while ((optc = getopt_long(argc, argv, GETOPTS, long_opts, NULL)) != -1) {
        switch(optc) {
            case 'h':
                print_help();
                return false;
                break;
            case 'i':
                set_config_param("input_file", strdup(optarg));
                break;
            case 'c':
                set_config_param("config_filepath", strdup(optarg));
                read_config();
                break;
            case 'd':
                set_config_param("input_file_directory", strdup(optarg));
                break;
            case 'p':
                set_config_param("papi_config_file", strdup(optarg));
                break;
            case 'o':
                set_config_param("output_file_prefix", strdup(optarg));
                break;
            case 'm':
                conf.mesh_duplicate_count = atoi(optarg);
                break;
            case 'g':
                conf.num_cycles = atoi(optarg);
                break;
            case 'v':
                conf.validate_result = true;
            case '\0':
                break;
            default:
                printf("Unknown command line parameter '%c'\n", optc);
        }
    }

    return true;
}

void print_config() {
    printf("---------- CONFIG ---------------------\n");
    printf("Config file: '%s'\n",          conf.config_filepath     ==NULL ? "NOT SET" : conf.config_filepath);
    printf("Input file: '%s'\n",           conf.input_file          ==NULL ? "NOT SET" : conf.input_file);
    printf("Input file directory: '%s'\n", conf.input_file_directory==NULL ? "NOT SET" : conf.input_file_directory);
    printf("PAPI config file: '%s'\n",     conf.papi_config_file    ==NULL ? "NOT SET" : conf.papi_config_file);
    printf("Output file prefix: '%s'\n",   conf.output_file_prefix  ==NULL ? "NOT SET" : conf.output_file_prefix);

    printf("OpenMP number of threads: %d\n", conf.omp_num_threads);

    printf("Mesh duplicate multiplier: %dx\n", conf.mesh_duplicate_count);
    printf("Multigrid V-cycles: %d\n", conf.num_cycles);

    if (conf.output_variables) {
        printf("Will write Euler equation variable values\n");
    }

    printf("---------------------------------------\n");
}

void print_help(void)
{
    fprintf(stderr, "MG-CFD instructions\n\n");
    fprintf(stderr, "Usage: euler3d_cpu.b [OPTIONS] \n\n");
    fprintf(stderr, "  -h, --help                       Print help\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "CRITICAL ARGUMENTS\n");
    fprintf(stderr, "  One of these must be set:\n");
    fprintf(stderr, "  -i, --input-file=FILEPATH        Multigrid input grid (.dat file)\n");
    fprintf(stderr, "  -c, --config-filepath=FILEPATH   Config file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONAL ARGUMENTS\n");
    fprintf(stderr, "  -d, --input-directory=DIRPATH    Directory path to input files\n");
    fprintf(stderr, "  -o, --output-file-prefix=STRING  String to prepend to output filenames\n");
    fprintf(stderr, "  -i, --papi-config-file=FILEPATH  PAPI events to monitor\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -m, --mesh-duplicate-count=INT   Number of times to duplicate mesh\n");
    fprintf(stderr, "  -g, --num-cycles=INT             Number of multigrid V-cycles\n");
    fprintf(stderr, "  -v, --validate-result            Check final state against pre-calculated solution\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "DEBUGGING ARGUMENTS\n");
    fprintf(stderr, "  --output-variables               Write Euler equation variable values to file\n");
    fprintf(stderr, "  --output-fluxes                  Write flux accumulations to file\n");
    fprintf(stderr, "  --output-step-factors            Write step factors to file\n");
}
