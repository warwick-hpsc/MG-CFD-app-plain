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
  #undef COLOURED_CONFLICT_AVOIDANCE
#endif

#endif