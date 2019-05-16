#ifndef COMMON_H
#define COMMON_H

#include <fstream>
#include <math.h>
#include <sstream>
#include <stdarg.h>
#include <stdlib.h>
#include <string>

#if defined OMP
#include <omp.h>
#endif

#include "config.h"
#include "const.h"
#include "definitions.h"
#include "globals.h"

std::string& ltrim(std::string &s);
std::string& rtrim(std::string &s);
std::string& trim(std::string &s);

const double smoothing_coefficient = double(0.2f);

const int point_fields[NVAR] = { VAR_DENSITY, VAR_MOMENTUMX, VAR_MOMENTUMY, VAR_MOMENTUMZ, VAR_DENSITY_ENERGY };

inline void log(const char *format, ...)
{
    #ifdef LOG
        va_list argp;
        #if defined OMP
            printf("Thread %d: ", omp_get_thread_num());
        #endif

        va_start(argp, format);
        vprintf(format, argp);
        va_end(argp);
        printf("\n");
    #endif
}

inline void zero_array(int nelr, double* array)
{
    for(int i=0; i<nelr; i++)
    {
        array[i] = 0.0;
    }
}

inline void zero_edges(
    int first_edge,
    int nedges,
    edge* restrict edge_variables)
{
    #ifdef OMP
        #pragma omp parallel for
    #endif
    for(int e=first_edge; e<(first_edge+nedges); e++)
    {
        for (int v=0; v<NVAR; v++) {
            edge_variables[e*NVAR + v].a = double(0.0);
            edge_variables[e*NVAR + v].b = double(0.0);
        }
    }
}

template <typename T>
std::string number_to_string(T number)
{
    std::ostringstream ss;
    ss << number;
    return ss.str();
}

template<typename T>
T* alloc(int N)
{
    #ifdef __ICC
        return new T[N];
    #else
        return (T*)malloc(N*sizeof(T));
    #endif
}

template<typename T>
void dealloc(T *restrict array)
{
    #ifdef __ICC
        delete[] array;
    #else
        free(array);
    #endif
}

template <typename T>
inline void copy(T *restrict dst, T *restrict src, int N)
{
    log("copy()");

    #ifdef OMP
        #pragma omp parallel for
    #endif
    for(int i=0; i<N; i++)
    {
        dst[i] = src[i];
    }
}

inline std::string get_cpu_model_name()
{
    std::ifstream file("/proc/cpuinfo");
    std::string file_line;

    if (!file.is_open()) {
        // fprintf(stderr, "ERROR: Could not open '/proc/cpuinfo'\n");
        // exit(EXIT_FAILURE);
        printf("WARNING: Could not open '/proc/cpuinfo'\n");
    }
    else {
        while (std::getline(file, file_line))
        {
            if (file_line.compare(0, 10, "model name") == 0)
            {
                std::size_t start = file_line.find(":");
                if (start == std::string::npos) {
                    // fprintf(stderr, "ERROR: Could not find ':' in model name line\n");
                    // exit(EXIT_FAILURE);
                    printf("WARNING: Could not find ':' in model name line\n");
                    break;
                }
                std::string model_name = file_line.substr(start+1, file_line.length());
                model_name = ltrim(model_name);
                return model_name;
            }
        }
    }
    return std::string("");
}

inline int compare_two_ints(const void * a, const void * b) {
    return ( *(int*)a - *(int*)b );
}

inline int compare_two_edges(edge_neighbour e1, edge_neighbour e2) {
    if (e1.a < e2.a) return true;
    if (e1.a > e2.a) return false;
    if (e1.b < e2.b) return true;
    if (e1.b > e2.b) return false;
    if (e1.x < e2.x) return true;
    if (e1.x > e2.x) return false;
    if (e1.y < e2.y) return true;
    if (e1.y > e2.y) return false;
    if (e1.z < e2.z) return true;
    if (e1.z > e2.z) return false;
    return false;
}

inline bool validate_edge_array(
    edge* restrict edge_variables,
    int length)
{
    for (int i=0; i<length; i++) {
        for (int v=0; v<NVAR; v++) {
            const double val_a = edge_variables[i*NVAR + v].a;
            const double val_b = edge_variables[i*NVAR + v].b;
            if (isnan(val_a) || isnan(-val_a)) {
                fprintf(stderr, "NaN detected at edge_variables[i=%d][v=%d].a\n", i, v);
                return false;
            }
            if (isinf(val_a)) {
                fprintf(stderr, "Infinity detected at edge_variables[i=%d][v=%d].a\n", i, v);
                return false;
            }
            if (isnan(val_b) || isnan(-val_b)) {
                fprintf(stderr, "NaN detected at edge_variables[i=%d][v=%d].b\n", i, v);
                return false;
            }
            if (isinf(val_b)) {
                fprintf(stderr, "Infinity detected at edge_variables[i=%d][v=%d].b\n", i, v);
                return false;
            }
        }
    }
    return true;
}

template <typename T>
void permute_array_range(T *restrict array, int *restrict mapping, int map_offset, int map_size)
{
    T* permutation = alloc<T>(map_size);

    bool* free_slots = alloc<bool>(map_size);
    for (int i=0; i<map_size; i++) free_slots[i] = true;

    for (int m=0; m<map_size; m++) {
        int from_idx = m + map_offset;
        int to_idx = mapping[m];
        permutation[to_idx] = array[from_idx];
        free_slots[to_idx] = false;
    }

    // // Validate mapping:
    // for (int i=0; i<map_size; i++) {
    //     if (free_slots[i]) {
    //         fprintf(stderr, "ERROR: No mapping made to destination %d\n", i+map_offset);
    //         DEBUGGABLE_ABORT
    //     }
    // }

    for (int i=0; i<map_size; i++) {
        array[map_offset+i] = permutation[i];
    }

    dealloc<T>(permutation);
    dealloc<bool>(free_slots);
}

bool file_exists(const char* filepath);

#ifdef OMP
inline void openmp_distribute_loop_iterations(int* start, int* end) {
    const int end_copy = *end;
    const int tid = omp_get_thread_num();
    const int num_t = omp_get_num_threads();
    const int thread_range = ((*end)-(*start)+num_t-1)/num_t;

    *start = (*start) + thread_range*tid;
    *end = (*start)+thread_range;
    if ((*end) > end_copy) *end = end_copy;
}
#endif

#endif