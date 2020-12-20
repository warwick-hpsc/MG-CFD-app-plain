#ifndef LOOP_STATS_H
#define LOOP_STATS_H

#include "common.h"

extern int iters_monitoring_state;

void init_iters();

void record_iters(long loop_start, long loop_end);

void dump_loop_stats_to_file();

#endif 
