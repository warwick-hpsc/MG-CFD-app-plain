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

#if defined _CRAYC
    #define UNROLL_LOOP(X) _Pragma( MACRO_TO_STR( _CRI unroll X ) )
    #define NOUNROLL _Pragma( MACRO_TO_STR( _CRI nounroll ) )
#elif defined __GNUC__ && ! defined __clang__
    #define UNROLL_LOOP(X) _Pragma( MACRO_TO_STR( GCC unroll X ) )
    #define NOUNROLL _Pragma( MACRO_TO_STR( GCC unroll 1 ) )
#else
    #define UNROLL_LOOP(X) _Pragma( MACRO_TO_STR( unroll X ) )
    #define NOUNROLL _Pragma( MACRO_TO_STR( unroll 1 ) )
#endif

#ifdef _CRAYC
    // Cray generates incorrect code if OMP SIMD pragma used:
    // I should submit a bug report ...
    #define NOSIMD _Pragma( MACRO_TO_STR( _CRI novector ) )
#elif defined __clang__
    #define NOSIMD _Pragma( MACRO_TO_STR( clang loop vectorize(disable) ) )
#else
    #define NOSIMD _Pragma( MACRO_TO_STR( omp simd safelen(1) ) )
#endif

// #ifdef __clang__
#if defined __clang__ || defined __arm__
    // With Clang, need to explicitly disable unrolling in case it interferes with SIMD:
    #define SIMD_LOOP(X)   _Pragma( MACRO_TO_STR( clang loop vectorize_width(X) interleave(disable) ) )
    #define SIMD_LOOP_AUTO _Pragma( MACRO_TO_STR( clang loop vectorize(enable)  interleave(disable) ) )
#else
    #define SIMD_LOOP(X)   _Pragma( MACRO_TO_STR( omp simd simdlen(DBLS_PER_SIMD) ) )
    #define SIMD_LOOP_AUTO _Pragma( MACRO_TO_STR( omp simd ) )
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
