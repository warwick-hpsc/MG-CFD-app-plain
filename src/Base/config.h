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

#ifndef CONFIG_H
#define CONFIG_H

#include <getopt.h>

typedef struct {
	char* config_filepath;

	char* input_file;
	char* input_file_directory;
	char* papi_config_file;
	char* output_file_prefix;

	int mesh_duplicate_count;
	int num_cycles;
	int omp_num_threads;

	bool renumber_mesh;

	bool validate_result;

	bool measure_mem_bound;
	bool perform_uns_compute;

	bool output_variables;
	bool output_old_variables;
	bool output_step_factors;
	bool output_edge_fluxes;
	bool output_fluxes;
	bool output_volumes;
} config;

extern config conf;

void set_config_defaults();

void read_config();

void print_config();

void print_help();

bool parse_arguments(int argc, char** argv);

#endif
