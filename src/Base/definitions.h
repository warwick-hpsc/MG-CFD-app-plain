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

#if defined SIMD && !defined DBLS_PER_SIMD
    #error "DBLS_PER_SIMD not defined, necessary for SIMD"
#endif

#define INTS_PER_VECTOR (DBLS_PER_SIMD*8)

// Number of doubles in a cache-line
#define DBLS_IN_CLINE (64/sizeof(double))

// CPP language does not have C's "restrict" keyword, but it has
// "__restrict__" which does the same thing:
#ifdef __ICC
#define restrict __restrict__
#elif defined __GNUC__
#define restrict __restrict__
#elif defined __clang__
#define restrict __restrict__
#endif

#ifdef __ICC
#define __assume_int_perfect_mult(x,m) __assume((x)%(m)==0);
#elif defined __GNUC__
#define __assume_int_perfect_mult(x,m) x = (x) & ~((m)-1);
#else
#define __assume_int_perfect_mult(x,m)
#endif

#define IS_ALIGNED(PTR, BC) \
    (((uintptr_t)(const void *)(PTR)) % (BC) == 0)

#ifdef __ICC
    #define DECLARE_INT_ALIGNED(X) __assume(X%DBLS_IN_CLINE==0)
#else
    #define DECLARE_INT_ALIGNED(X) 
#endif

#ifdef __ICC
    #define DECLARE_PTR_ALIGNED(X) __assume_aligned(X, 64)
#else
    #define DECLARE_PTR_ALIGNED(X)
#endif

#ifdef __ICC
    #define IVDEP ivdep
    // #define IVDEP vector always
#else
    #define IVDEP GCC ivdep
    // #define IVDEP vector always
#endif

#define DEBUGGABLE_ABORT fprintf(stderr, "%s:%d\n", __FILE__, __LINE__); fflush(stderr); fflush(stdout); exit(EXIT_FAILURE);

//////////////////////
//////////////////////
//////////////////////

namespace Filetype
{
    enum Filetype {
        Unknown,
        HDF5,
        Old_style,
    };
}

struct double3 { double x, y, z; };
struct edge { double a, b; };
struct edge_neighbour { int a, b; double x, y, z; };

typedef struct {
    int length;
    int* data;
} ivector;

typedef struct {
    unsigned long size,resident,share,text,lib,data,dt;
} statm_t;

#endif
