#ifndef PAPI_FUNCS_H
#define PAPI_FUNCS_H

#ifdef PAPI

#include <papi.h>

#include "common.h"

inline unsigned long omp_get_thread_num_ul() {
    #if defined OMP
        return (unsigned long)omp_get_thread_num();
    #else
        return 0;
    #endif
}

void start_papi();

void stop_papi();

void init_papi();

void load_papi_events();

void dump_papi_counters_to_file();

#endif

#endif
