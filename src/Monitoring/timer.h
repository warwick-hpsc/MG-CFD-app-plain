#ifndef TIMER_H
#define TIMER_H

#ifdef TIME

#include <sys/time.h>

#include "common.h"

void init_timers();

void start_timer();

void stop_timer();

void dump_timers_to_file();

#endif

#endif 
