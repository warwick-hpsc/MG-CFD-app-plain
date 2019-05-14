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

#include <stdlib.h>
#include <cmath>

#include "mg.h"
#include "papi_funcs.h"
#include "timer.h"
#include "loop_stats.h"

void down_residuals(
    double *restrict residuals1, 
    int nel1, 
    double *restrict variables2, 
    double *restrict residuals2, 
    int nel2, 
    int *restrict mapping, 
    int mgc, 
    double3 *restrict coords1, 
    double3 *restrict coords2)
{
    log("down()");
    current_kernel = DOWN;

    int loop_start = 0;
    int loop_end = mgc;

    #ifdef OMP
        #pragma omp parallel firstprivate(loop_start, loop_end)
        {
            openmp_distribute_loop_iterations(&loop_start, &loop_end);
    #endif

    #ifdef PAPI
    start_papi();
    #endif
    #ifdef TIME
    start_timer();
    #endif
    for(int i=loop_start; i<loop_end; i++)
    {
        const int p1 = mapping[i];
        
        //1. Calculate dx, dy, dz, dm
        double dx = fabs(coords2[i].x - coords1[p1].x);
        double dy = fabs(coords2[i].y - coords1[p1].y);
        double dz = fabs(coords2[i].z - coords1[p1].z);
        double dm = sqrt(dx*dx + dy*dy + dz*dz);

        const int p_idx2  = NVAR*i + VAR_DENSITY;
        const int mx_idx2 = NVAR*i + VAR_MOMENTUMX;
        const int my_idx2 = NVAR*i + VAR_MOMENTUMY;
        const int mz_idx2 = NVAR*i + VAR_MOMENTUMZ;
        const int pe_idx2 = NVAR*i + VAR_DENSITY_ENERGY;

        const int p_idx1  = NVAR*p1 + VAR_DENSITY;
        const int mx_idx1 = NVAR*p1 + VAR_MOMENTUMX;
        const int my_idx1 = NVAR*p1 + VAR_MOMENTUMY;
        const int mz_idx1 = NVAR*p1 + VAR_MOMENTUMZ;
        const int pe_idx1 = NVAR*p1 + VAR_DENSITY_ENERGY;

        variables2[p_idx2]  -= (residuals1[p_idx1 ] - residuals2[p_idx2 ])*dm;
        variables2[mx_idx2] -= (residuals1[mx_idx1] - residuals2[mx_idx2])*dx;
        variables2[my_idx2] -= (residuals1[my_idx1] - residuals2[my_idx2])*dy;
        variables2[mz_idx2] -= (residuals1[mz_idx1] - residuals2[mz_idx2])*dz;
        variables2[pe_idx2] -= (residuals1[pe_idx1] - residuals2[pe_idx2])*dm;
    }
    #ifdef TIME
    stop_timer();
    #endif
    #ifdef PAPI
    stop_papi();
    #endif
    record_iters(loop_start, loop_end);

    #ifdef OMP
    }
    #endif
}

void up(
    double *restrict variables1, 
    double *restrict variables2, 
    int nel2, 
    int *restrict mapping, 
    int *restrict up_scratch, 
    int mgc)
{
    log("up()");
    current_kernel = UP;

    int loop_start, loop_end;

    // Zero the variables array:
    {
        // ... but only for nodes connected to level below.
        loop_start = 0;
        loop_end = mgc;

        #ifdef OMP
            #pragma omp parallel firstprivate(loop_start, loop_end)
            {
                openmp_distribute_loop_iterations(&loop_start, &loop_end);
        #endif

        #ifdef PAPI
        start_papi();
        #endif
        #ifdef TIME
        start_timer();
        #endif
        for(int i=loop_start; i<loop_end; i++)
        {
            int p2 = mapping[i];

            const int p_idx2  = NVAR*p2 + VAR_DENSITY;
            const int mx_idx2 = NVAR*p2 + VAR_MOMENTUMX;
            const int my_idx2 = NVAR*p2 + VAR_MOMENTUMY;
            const int mz_idx2 = NVAR*p2 + VAR_MOMENTUMZ;
            const int pe_idx2 = NVAR*p2 + VAR_DENSITY_ENERGY;

            variables2[p_idx2]  = 0.0;
            variables2[mx_idx2] = 0.0;
            variables2[my_idx2] = 0.0;
            variables2[mz_idx2] = 0.0;
            variables2[pe_idx2] = 0.0;
        }
        #ifdef TIME
        stop_timer();
        #endif
        #ifdef PAPI
        stop_papi();
        #endif

        #ifdef OMP
        }
        #endif
    }

    // Zero the up_scratch array:
    {
        #pragma omp parallel for
        for (int i=0; i<nel2; i++)
        {
            up_scratch[i] = 0;
        }
    }

    // Accumulate from level below:
    {
        #ifdef PAPI
        start_papi();
        #endif
        #ifdef TIME
        start_timer();
        #endif
        for(int i=0; i<mgc; i++)
        {
            int p2 = mapping[i];

            const int p_idx2  = NVAR*p2 + VAR_DENSITY;
            const int mx_idx2 = NVAR*p2 + VAR_MOMENTUMX;
            const int my_idx2 = NVAR*p2 + VAR_MOMENTUMY;
            const int mz_idx2 = NVAR*p2 + VAR_MOMENTUMZ;
            const int pe_idx2 = NVAR*p2 + VAR_DENSITY_ENERGY;

            const int p_idx1  = NVAR*i + VAR_DENSITY;
            const int mx_idx1 = NVAR*i + VAR_MOMENTUMX;
            const int my_idx1 = NVAR*i + VAR_MOMENTUMY;
            const int mz_idx1 = NVAR*i + VAR_MOMENTUMZ;
            const int pe_idx1 = NVAR*i + VAR_DENSITY_ENERGY;

            variables2[p_idx2]  += variables1[p_idx1];
            variables2[mx_idx2] += variables1[mx_idx1];
            variables2[my_idx2] += variables1[my_idx1];
            variables2[mz_idx2] += variables1[mz_idx1];
            variables2[pe_idx2] += variables1[pe_idx1];

            up_scratch[p2]++;
        }
        #ifdef TIME
        stop_timer();
        #endif
        #ifdef PAPI
        stop_papi();
        #endif
    }
    
    // Apply average to accumulation:
    {
        loop_start = 0;
        loop_end = nel2;
        #ifdef OMP
            #pragma omp parallel firstprivate(loop_start, loop_end)
            {
                openmp_distribute_loop_iterations(&loop_start, &loop_end);
        #endif

        #ifdef PAPI
        start_papi();
        #endif
        #ifdef TIME
        start_timer();
        #endif
        for(int i=loop_start; i<loop_end; i++)
        {
            double average = up_scratch[i]==0 ? 1.0 : 1.0 / (double)up_scratch[i];

            const int p_idx2  = NVAR*i + VAR_DENSITY;
            const int mx_idx2 = NVAR*i + VAR_MOMENTUMX;
            const int my_idx2 = NVAR*i + VAR_MOMENTUMY;
            const int mz_idx2 = NVAR*i + VAR_MOMENTUMZ;
            const int pe_idx2 = NVAR*i + VAR_DENSITY_ENERGY;

            variables2[p_idx2]  *= average;
            variables2[mx_idx2] *= average;
            variables2[my_idx2] *= average;
            variables2[mz_idx2] *= average;
            variables2[pe_idx2] *= average;
        }
        #ifdef TIME
        stop_timer();
        #endif
        #ifdef PAPI
        stop_papi();
        #endif
        record_iters(loop_start, loop_end);

        #ifdef OMP
        }
        #endif
    }
}
