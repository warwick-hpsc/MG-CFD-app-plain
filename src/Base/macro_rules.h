#ifndef MACRO_RULES
#define MACRO_RULES

#define VALUE_TO_STRING(x) #x
#define VALUE(x) VALUE_TO_STRING(x)
#define VAR_NAME_VALUE(var) #var "=" VALUE( (var) )
#define XMACRO_TO_STR(s) MACRO_TO_STR( s )
#define MACRO_TO_STR(s) #s
#define DO_PRAGMA(x) _Pragma ( #x )

#if defined BIN_COLOURED_VECTORS || defined BIN_COLOURED_CONTIGUOUS
  #define COLOURED_CONFLICT_AVOIDANCE 1
#elif defined COLOURED_CONFLICT_AVOIDANCE
  #pragma message("Colouring scheme not specified, disabling coloured CA")
  #undef COLOURED_CONFLICT_AVOIDANCE
#endif

#if defined MANUAL_GATHER
	#ifndef MANUAL_SCATTER
		#define MANUAL_SCATTER
	#endif
#endif

#if defined MANUAL_GATHER || defined MANUAL_SCATTER
  #ifndef MANUAL_WIDTH
    // #define MANUAL_WIDTH 512
    #ifdef __clang__
      // Despite attempts to use aligned arrays, Clang-generated binary is 
      // determined to perform 2 serial loop iterations for each SIMD-loop full pass. 
      // So to get any SIMD execution, need at least 2+simd_width. 
      // Then Amadhl's law applies, as those serial iterations prevent 
      // full 4x SIMD speedup. With 3x SIMD iterations, get 2.8x fewer 
      // iterations overall, theoretically a sweet spot (as going wider 
      // potentially hits  memory limits with all that batched packing/unpacking):
      #define MANUAL_WIDTH ((DBLS_PER_SIMD*3)+2)
    #else
      #define MANUAL_WIDTH DBLS_PER_SIMD
    #endif
  #endif
#endif

#if defined SIMD && !(defined FLUX_FISSION)
    #ifdef COLOURED_CONFLICT_AVOIDANCE
    #elif defined MANUAL_SCATTER
    #elif defined __AVX512CD__ && defined __ICC
		#pragma message("Enabling USE_AVX512CD flag")
        #define USE_AVX512CD 1
    #endif
#endif

#endif
