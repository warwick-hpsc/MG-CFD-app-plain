#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include "const.h"

//////////////////////
/// Macros:
//////////////////////


#define VALUE_TO_STRING(x) #x
#define VALUE(x) VALUE_TO_STRING(x)
#define VAR_NAME_VALUE(var) #var "=" VALUE( (var) )
#define XMACRO_TO_STR(s) MACRO_TO_STR( s )
#define MACRO_TO_STR(s) #s
#define DO_PRAGMA(x) _Pragma ( #x )


#ifndef DBLS_PER_SIMD
    #ifndef SIMD
        #define DBLS_PER_SIMD 1
    #endif
#endif

// CPP language does not have C's "restrict" keyword, but it has
// "__restrict__" which does the same thing:
#ifdef __ICC
#define restrict __restrict__
#elif defined __GNUC__
#define restrict __restrict__
#elif defined __clang__
#define restrict __restrict__
#endif

#if defined __clang__ || defined __GNUC__
#define FORCE_INLINE __attribute__((always_inline))
#else
#define FORCE_INLINE
#endif

#if defined __GNUC__ && ! defined __clang__
#define UNROLL_LOOP_FULLY _Pragma("GCC unroll 1000")
#else
#define UNROLL_LOOP_FULLY _Pragma("unroll")
#endif

#define DEBUGGABLE_ABORT fprintf(stderr, "%s:%d\n", __FILE__, __LINE__); fflush(stderr); fflush(stdout); exit(EXIT_FAILURE);

//////////////////////
//////////////////////
//////////////////////

struct double3 { double x, y, z; };
struct edge { double a, b; };
struct edge_neighbour { long a, b; double x, y, z; };

typedef struct {
    int length;
    int* data;
} ivector;

typedef struct {
    unsigned long size,resident,share,text,lib,data,dt;
} statm_t;

#endif
