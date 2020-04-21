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

#include "mg_loops.h"

#include "papi_funcs.h"
#include "timer.h"
#include "loop_stats.h"

void mg_restrict(
    double *restrict variables1, 
    double *restrict variables2, 
    long nel2, 
    long *restrict mapping, 
    long *restrict up_scratch, 
    long mgc)
{
    log("restrict()");
    current_kernel = RESTRICT;

    long loop_start, loop_end;

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
        record_iters(loop_start, loop_end);

        for(long i=loop_start; i<loop_end; i++)
        {
            long p2 = mapping[i];

            const long p_idx2  = NVAR*p2 + VAR_DENSITY;
            const long mx_idx2 = NVAR*p2 + VAR_MOMENTUMX;
            const long my_idx2 = NVAR*p2 + VAR_MOMENTUMY;
            const long mz_idx2 = NVAR*p2 + VAR_MOMENTUMZ;
            const long pe_idx2 = NVAR*p2 + VAR_DENSITY_ENERGY;

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
        for (long i=0; i<nel2; i++)
        {
            up_scratch[i] = 0;
        }
    }

    // Accumulate from level below:
    {
        loop_start = 0;
        loop_end = mgc;
        #if defined OMP && defined OMP_SCATTERS
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
        record_iters(loop_start, loop_end);

        for(long i=loop_start; i<loop_end; i++)
        {
            long p2 = mapping[i];

            const long p_idx2  = NVAR*p2 + VAR_DENSITY;
            const long mx_idx2 = NVAR*p2 + VAR_MOMENTUMX;
            const long my_idx2 = NVAR*p2 + VAR_MOMENTUMY;
            const long mz_idx2 = NVAR*p2 + VAR_MOMENTUMZ;
            const long pe_idx2 = NVAR*p2 + VAR_DENSITY_ENERGY;

            const long p_idx1  = NVAR*i + VAR_DENSITY;
            const long mx_idx1 = NVAR*i + VAR_MOMENTUMX;
            const long my_idx1 = NVAR*i + VAR_MOMENTUMY;
            const long mz_idx1 = NVAR*i + VAR_MOMENTUMZ;
            const long pe_idx1 = NVAR*i + VAR_DENSITY_ENERGY;

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

        #if defined OMP && defined OMP_SCATTERS
        }
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
        record_iters(loop_start, loop_end);

        for(long i=loop_start; i<loop_end; i++)
        {
            double average = up_scratch[i]==0 ? 1.0 : 1.0 / (double)up_scratch[i];

            const long p_idx2  = NVAR*i + VAR_DENSITY;
            const long mx_idx2 = NVAR*i + VAR_MOMENTUMX;
            const long my_idx2 = NVAR*i + VAR_MOMENTUMY;
            const long mz_idx2 = NVAR*i + VAR_MOMENTUMZ;
            const long pe_idx2 = NVAR*i + VAR_DENSITY_ENERGY;

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

        #ifdef OMP
        }
        #endif
    }
}

void prolong(
    double *restrict variables1, 
    double *restrict variables2, 
    long *restrict mapping, 
    long mgc, 
    double3 *restrict coords1, 
    double3 *restrict coords2)
{
    // This is the original 'prolong' operator added by my predecessor. 
    // I think it is mathematically flawed, quickly corrupting the solution.
    // Attempts to fix are made in later functions.

    log("prolong()");
    current_kernel = PROLONG;

    long loop_start = 0;
    long loop_end = mgc;

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
    record_iters(loop_start, loop_end);

    for(long i=loop_start; i<loop_end; i++)
    {
        const long p1 = mapping[i];
        
        //1. Calculate dx, dy, dz, dm
        double dx = fabs(coords2[i].x - coords1[p1].x);
        double dy = fabs(coords2[i].y - coords1[p1].y);
        double dz = fabs(coords2[i].z - coords1[p1].z);
        double dm = sqrt(dx*dx + dy*dy + dz*dz);

        const long p_idx2  = NVAR*i + VAR_DENSITY;
        const long mx_idx2 = NVAR*i + VAR_MOMENTUMX;
        const long my_idx2 = NVAR*i + VAR_MOMENTUMY;
        const long mz_idx2 = NVAR*i + VAR_MOMENTUMZ;
        const long pe_idx2 = NVAR*i + VAR_DENSITY_ENERGY;

        const long p_idx1  = NVAR*p1 + VAR_DENSITY;
        const long mx_idx1 = NVAR*p1 + VAR_MOMENTUMX;
        const long my_idx1 = NVAR*p1 + VAR_MOMENTUMY;
        const long mz_idx1 = NVAR*p1 + VAR_MOMENTUMZ;
        const long pe_idx1 = NVAR*p1 + VAR_DENSITY_ENERGY;

        variables2[p_idx2]  -= (variables1[p_idx1 ] - variables2[p_idx2 ])*dm;
        variables2[mx_idx2] -= (variables1[mx_idx1] - variables2[mx_idx2])*dx;
        variables2[my_idx2] -= (variables1[my_idx1] - variables2[my_idx2])*dy;
        variables2[mz_idx2] -= (variables1[mz_idx1] - variables2[mz_idx2])*dz;
        variables2[pe_idx2] -= (variables1[pe_idx1] - variables2[pe_idx2])*dm;
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

void prolong_residuals(
    double *restrict residuals1, 
    // Depending on MG configuration variables2 and residuals2
    // may point to the same array, so cannot use 'restrict' 
    // qualifier:
    double *variables2, 
    double *residuals2, 
    long *restrict mapping, 
    long mgc)
{
    log("prolong_residuals()");
    current_kernel = PROLONG;

    long loop_start = 0;
    long loop_end = mgc;

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
    record_iters(loop_start, loop_end);

    for(long i=loop_start; i<loop_end; i++)
    {
        const long p1 = mapping[i];
        
        const long p_idx2  = NVAR*i + VAR_DENSITY;
        const long mx_idx2 = NVAR*i + VAR_MOMENTUMX;
        const long my_idx2 = NVAR*i + VAR_MOMENTUMY;
        const long mz_idx2 = NVAR*i + VAR_MOMENTUMZ;
        const long pe_idx2 = NVAR*i + VAR_DENSITY_ENERGY;

        const long p_idx1  = NVAR*p1 + VAR_DENSITY;
        const long mx_idx1 = NVAR*p1 + VAR_MOMENTUMX;
        const long my_idx1 = NVAR*p1 + VAR_MOMENTUMY;
        const long mz_idx1 = NVAR*p1 + VAR_MOMENTUMZ;
        const long pe_idx1 = NVAR*p1 + VAR_DENSITY_ENERGY;

        variables2[p_idx2]  += (residuals2[p_idx1 ] - residuals1[p_idx2 ]);
        variables2[mx_idx2] += (residuals2[mx_idx1] - residuals1[mx_idx2]);
        variables2[my_idx2] += (residuals2[my_idx1] - residuals1[my_idx2]);
        variables2[mz_idx2] += (residuals2[mz_idx1] - residuals1[mz_idx2]);
        variables2[pe_idx2] += (residuals2[pe_idx1] - residuals1[pe_idx2]);
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

void prolong_interpolate(
    double *restrict variables1, 
    long nel1, 
    double *restrict variables2, 
    long *restrict mapping, 
    long mgc, 
    double3 *restrict coords1, 
    double3 *restrict coords2)
{
    log("prolong_interpolate()");
    current_kernel = PROLONG;

    long loop_start = 0;
    long loop_end = mgc;

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
    record_iters(loop_start, loop_end);

    for(long i=loop_start; i<loop_end; i++)
    {
        const long p2 = i;

        const long p1  = mapping[p2];

        // It would be better to get actual neighbours 
        // of node p1, but for now am just experimenting:
        long p1a, p1b;
        if (p1==0) {
            p1a = p1+1;
            p1b = p1+2;
        }
        else if (p1==(nel1-1)) {
            p1a = p1-2;
            p1b = p1-1;
        }
        else {
            p1a = p1-1;
            p1b = p1+1;
        }

        double dx = fabs(coords2[p2].x - coords1[p1].x);
        double dy = fabs(coords2[p2].y - coords1[p1].y);
        double dz = fabs(coords2[p2].z - coords1[p1].z);
        double dm = sqrt(dx*dx + dy*dy + dz*dz);

        double dx_a = fabs(coords2[p2].x - coords1[p1a].x);
        double dy_a = fabs(coords2[p2].y - coords1[p1a].y);
        double dz_a = fabs(coords2[p2].z - coords1[p1a].z);
        double dm_a = sqrt(dx_a*dx_a + dy_a*dy_a + dz_a*dz_a);

        double dx_b = fabs(coords2[p2].x - coords1[p1b].x);
        double dy_b = fabs(coords2[p2].y - coords1[p1b].y);
        double dz_b = fabs(coords2[p2].z - coords1[p1b].z);
        double dm_b = sqrt(dx_b*dx_b + dy_b*dy_b + dz_b*dz_b);

        // double dm_sum = dm + dm_a + dm_b;

        double p2_factor, p2a_factor, p2b_factor, w_sum;
        if (dm == 0.0) {
            // Simply transfer the value down:
            p2_factor = 1.0;
            p2a_factor = 0.0;
            p2b_factor = 0.0;
            w_sum = 1.0;
        }
        else {
            // Interpolate:
            // p2_factor  = (dm  ==0.0) ? 0.0 : dm  /dm_sum;
            // p2a_factor = (dm_a==0.0) ? 0.0 : dm_a/dm_sum;
            // p2b_factor = (dm_b==0.0) ? 0.0 : dm_b/dm_sum;
            p2_factor  = (dm  ==0.0) ? 0.0 : 1.0 / dm;
            p2a_factor = (dm_a==0.0) ? 0.0 : 1.0 / dm_a;
            p2b_factor = (dm_b==0.0) ? 0.0 : 1.0 / dm_b;
            w_sum = p2_factor + p2a_factor + p2b_factor;
        }

        const long p_idx2  = NVAR*p2 + VAR_DENSITY;
        const long mx_idx2 = NVAR*p2 + VAR_MOMENTUMX;
        const long my_idx2 = NVAR*p2 + VAR_MOMENTUMY;
        const long mz_idx2 = NVAR*p2 + VAR_MOMENTUMZ;
        const long pe_idx2 = NVAR*p2 + VAR_DENSITY_ENERGY;

        const long p_idx1  = NVAR*p1 + VAR_DENSITY;
        const long mx_idx1 = NVAR*p1 + VAR_MOMENTUMX;
        const long my_idx1 = NVAR*p1 + VAR_MOMENTUMY;
        const long mz_idx1 = NVAR*p1 + VAR_MOMENTUMZ;
        const long pe_idx1 = NVAR*p1 + VAR_DENSITY_ENERGY;

        const long p_idx1a  = NVAR*p1a + VAR_DENSITY;
        const long mx_idx1a = NVAR*p1a + VAR_MOMENTUMX;
        const long my_idx1a = NVAR*p1a + VAR_MOMENTUMY;
        const long mz_idx1a = NVAR*p1a + VAR_MOMENTUMZ;
        const long pe_idx1a = NVAR*p1a + VAR_DENSITY_ENERGY;

        const long p_idx1b  = NVAR*p1b + VAR_DENSITY;
        const long mx_idx1b = NVAR*p1b + VAR_MOMENTUMX;
        const long my_idx1b = NVAR*p1b + VAR_MOMENTUMY;
        const long mz_idx1b = NVAR*p1b + VAR_MOMENTUMZ;
        const long pe_idx1b = NVAR*p1b + VAR_DENSITY_ENERGY;

        variables2[p_idx2] = p2_factor*variables1[p_idx1] 
                           + p2a_factor*variables1[p_idx1a]
                           + p2b_factor*variables1[p_idx1b];

        variables2[mx_idx2] = p2_factor*variables1[mx_idx1] 
                            + p2a_factor*variables1[mx_idx1a]
                            + p2b_factor*variables1[mx_idx1b];
                           
        variables2[my_idx2] = p2_factor*variables1[my_idx1] 
                            + p2a_factor*variables1[my_idx1a]
                            + p2b_factor*variables1[my_idx1b];
                           
        variables2[mz_idx2] = p2_factor*variables1[mz_idx1] 
                            + p2a_factor*variables1[mz_idx1a]
                            + p2b_factor*variables1[mz_idx1b];

        variables2[pe_idx2] = p2_factor*variables1[pe_idx1] 
                            + p2a_factor*variables1[pe_idx1a]
                            + p2b_factor*variables1[pe_idx1b];

        variables2[p_idx2]  /= w_sum;
        variables2[mx_idx2] /= w_sum;
        variables2[my_idx2] /= w_sum;
        variables2[mz_idx2] /= w_sum;
        variables2[pe_idx2] /= w_sum;
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

void prolong_residuals_interpolate_crude(
    double *restrict residuals1, 
    long nel1, 
    double *restrict residuals2,
    double *restrict variables2, 
    long *restrict mapping, 
    long mgc, 
    double3 *restrict coords1, 
    double3 *restrict coords2)
{
    log("prolong_residuals_interpolate_crude()");
    current_kernel = PROLONG;

    long loop_start = 0;
    long loop_end = mgc;

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
    record_iters(loop_start, loop_end);

    for(long i=loop_start; i<loop_end; i++)
    {
        const long p2 = i;

        const long p1  = mapping[p2];

        // It would be better to get actual neighbours of node p1, but that 
        // requires generating a mapping of node -> N nearest neighbours 
        // (where N is constant to simplify mapping access, eg 3).
        long p1a, p1b;
        if (p1==0) {
            p1a = p1+1;
            p1b = p1+2;
        }
        else if (p1==(nel1-1)) {
            p1a = p1-2;
            p1b = p1-1;
        }
        else {
            p1a = p1-1;
            p1b = p1+1;
        }

        double dx = fabs(coords2[p2].x - coords1[p1].x);
        double dy = fabs(coords2[p2].y - coords1[p1].y);
        double dz = fabs(coords2[p2].z - coords1[p1].z);
        double dm = sqrt(dx*dx + dy*dy + dz*dz);

        double dx_a = fabs(coords2[p2].x - coords1[p1a].x);
        double dy_a = fabs(coords2[p2].y - coords1[p1a].y);
        double dz_a = fabs(coords2[p2].z - coords1[p1a].z);
        double dm_a = sqrt(dx_a*dx_a + dy_a*dy_a + dz_a*dz_a);

        double dx_b = fabs(coords2[p2].x - coords1[p1b].x);
        double dy_b = fabs(coords2[p2].y - coords1[p1b].y);
        double dz_b = fabs(coords2[p2].z - coords1[p1b].z);
        double dm_b = sqrt(dx_b*dx_b + dy_b*dy_b + dz_b*dz_b);

        // double dm_sum = dm + dm_a + dm_b;

        double p2_factor, p2a_factor, p2b_factor, w_sum;
        if (dm == 0.0) {
            // Simply transfer the value down:
            p2_factor = 1.0;
            p2a_factor = 0.0;
            p2b_factor = 0.0;
            w_sum = 1.0;
        }
        else {
            // Interpolate:
            // p2_factor  = (dm  ==0.0) ? 0.0 : dm  /dm_sum;
            // p2a_factor = (dm_a==0.0) ? 0.0 : dm_a/dm_sum;
            // p2b_factor = (dm_b==0.0) ? 0.0 : dm_b/dm_sum;
            p2_factor  = (dm  ==0.0) ? 0.0 : 1.0 / dm;
            p2a_factor = (dm_a==0.0) ? 0.0 : 1.0 / dm_a;
            p2b_factor = (dm_b==0.0) ? 0.0 : 1.0 / dm_b;
            w_sum = p2_factor + p2a_factor + p2b_factor;
        }

        const long p_idx2  = NVAR*p2 + VAR_DENSITY;
        const long mx_idx2 = NVAR*p2 + VAR_MOMENTUMX;
        const long my_idx2 = NVAR*p2 + VAR_MOMENTUMY;
        const long mz_idx2 = NVAR*p2 + VAR_MOMENTUMZ;
        const long pe_idx2 = NVAR*p2 + VAR_DENSITY_ENERGY;

        const long p_idx1  = NVAR*p1 + VAR_DENSITY;
        const long mx_idx1 = NVAR*p1 + VAR_MOMENTUMX;
        const long my_idx1 = NVAR*p1 + VAR_MOMENTUMY;
        const long mz_idx1 = NVAR*p1 + VAR_MOMENTUMZ;
        const long pe_idx1 = NVAR*p1 + VAR_DENSITY_ENERGY;

        const long p_idx1a  = NVAR*p1a + VAR_DENSITY;
        const long mx_idx1a = NVAR*p1a + VAR_MOMENTUMX;
        const long my_idx1a = NVAR*p1a + VAR_MOMENTUMY;
        const long mz_idx1a = NVAR*p1a + VAR_MOMENTUMZ;
        const long pe_idx1a = NVAR*p1a + VAR_DENSITY_ENERGY;

        const long p_idx1b  = NVAR*p1b + VAR_DENSITY;
        const long mx_idx1b = NVAR*p1b + VAR_MOMENTUMX;
        const long my_idx1b = NVAR*p1b + VAR_MOMENTUMY;
        const long mz_idx1b = NVAR*p1b + VAR_MOMENTUMZ;
        const long pe_idx1b = NVAR*p1b + VAR_DENSITY_ENERGY;


        // variables2[p_idx2]  += p2_factor *residuals1[p_idx1] 
        //                     +  p2a_factor*residuals1[p_idx1a]
        //                     +  p2b_factor*residuals1[p_idx1b];

        // variables2[mx_idx2] += p2_factor *residuals1[mx_idx1] 
        //                     +  p2a_factor*residuals1[mx_idx1a]
        //                     +  p2b_factor*residuals1[mx_idx1b];
                           
        // variables2[my_idx2] += p2_factor *residuals1[my_idx1] 
        //                     +  p2a_factor*residuals1[my_idx1a]
        //                     +  p2b_factor*residuals1[my_idx1b];
                           
        // variables2[mz_idx2] += p2_factor *residuals1[mz_idx1] 
        //                     +  p2a_factor*residuals1[mz_idx1a]
        //                     +  p2b_factor*residuals1[mz_idx1b];

        // variables2[pe_idx2] += p2_factor *residuals1[pe_idx1] 
        //                     +  p2a_factor*residuals1[pe_idx1a]
        //                     +  p2b_factor*residuals1[pe_idx1b];
        
        // variables2[p_idx2]  /= w_sum;
        // variables2[mx_idx2] /= w_sum;
        // variables2[my_idx2] /= w_sum;
        // variables2[mz_idx2] /= w_sum;
        // variables2[pe_idx2] /= w_sum;

        double E[NVAR];
        E[VAR_DENSITY]   = p2_factor *residuals1[p_idx1] 
                         + p2a_factor*residuals1[p_idx1a]
                         + p2b_factor*residuals1[p_idx1b];

        E[VAR_MOMENTUMX] = p2_factor *residuals1[mx_idx1] 
                         + p2a_factor*residuals1[mx_idx1a]
                         + p2b_factor*residuals1[mx_idx1b];
                           
        E[VAR_MOMENTUMY] = p2_factor *residuals1[my_idx1] 
                         + p2a_factor*residuals1[my_idx1a]
                         + p2b_factor*residuals1[my_idx1b];
                           
        E[VAR_MOMENTUMZ] = p2_factor *residuals1[mz_idx1] 
                         + p2a_factor*residuals1[mz_idx1a]
                         + p2b_factor*residuals1[mz_idx1b];

        E[VAR_DENSITY_ENERGY] = p2_factor *residuals1[pe_idx1] 
                              + p2a_factor*residuals1[pe_idx1a]
                              + p2b_factor*residuals1[pe_idx1b];

        E[VAR_DENSITY]  /= w_sum;
        E[VAR_MOMENTUMX] /= w_sum;
        E[VAR_MOMENTUMY] /= w_sum;
        E[VAR_MOMENTUMZ] /= w_sum;
        E[VAR_DENSITY_ENERGY] /= w_sum;
      
        variables2[p_idx2]  += residuals2[p_idx2]  - E[VAR_DENSITY];
        variables2[mx_idx2] += residuals2[mx_idx2] - E[VAR_MOMENTUMX];
        variables2[my_idx2] += residuals2[my_idx2] - E[VAR_MOMENTUMY];
        variables2[mz_idx2] += residuals2[mz_idx2] - E[VAR_MOMENTUMZ];
        variables2[pe_idx2] += residuals2[pe_idx2] - E[VAR_DENSITY_ENERGY];
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

void prolong_residuals_interpolate_proper(
    edge_neighbour *edges,
    long num_edges,
    double *restrict residuals1, 
    double *restrict residuals2,
    double *restrict variables2, 
    long nel2,
    long *restrict mapping, 
    double3 *restrict coords1, 
    double3 *restrict coords2)
{
    // Goal: interpolate residuals of level above to estimate residuals of level below.
    //
    // For each node that has the same coordinates as its MG node parent, 
    // the 'prolonged residual' is simply taken directly from the MG node. 
    //
    // For each other node N, the 'prolonged residuals' is the weighted average 
    // across N's MG node and MG nodes of N's neighbours, requiring an 
    // edge-based loop. The weight is 1.0/distance.

    log("prolong_residuals_interpolate_proper()");
    current_kernel = PROLONG;

    double* w_sums = alloc<double>(nel2);
    for (long i=0; i<nel2; i++) {
        w_sums[i] = 0.0;
    }

    double* res2_wavg = alloc<double>(nel2*NVAR);
    for (long i=0; i<nel2*NVAR; i++) {
        res2_wavg[i] = 0.0;
    }

    // a1 and b1 belong to level above (L+1); a2 and b2 belong to level below (L)

    // 1) Perform the summing stage of weighted average:
    long loop_start = 0;
    long loop_end = loop_start + num_edges;
    #if defined OMP && (defined OMP_SCATTERS)
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
    record_iters(loop_start, loop_end);

    for (long i=loop_start; i<loop_end; i++) {
        const long a2 = edges[i].a;
        const long a1 = mapping[a2];
        const double3 ca1 = coords1[a1];
        const double3 ca2 = coords2[a2];

        const long b2 = edges[i].b;
        const long b1 = mapping[b2];
        const double3 cb1 = coords1[b1];
        const double3 cb2 = coords2[b2];

        // Process a2:
        double dx_a1a2 = ca2.x - ca1.x;
        double dy_a1a2 = ca2.y - ca1.y;
        double dz_a1a2 = ca2.z - ca1.z;
        if (dx_a1a2 == 0.0 && dy_a1a2 == 0.0 && dz_a1a2 == 0.0) {
            // a2 == a1:
            res2_wavg[a2*NVAR + VAR_DENSITY]        = residuals1[a1*NVAR + VAR_DENSITY];
            res2_wavg[a2*NVAR + VAR_MOMENTUMX]      = residuals1[a1*NVAR + VAR_MOMENTUMX];
            res2_wavg[a2*NVAR + VAR_MOMENTUMY]      = residuals1[a1*NVAR + VAR_MOMENTUMY];
            res2_wavg[a2*NVAR + VAR_MOMENTUMZ]      = residuals1[a1*NVAR + VAR_MOMENTUMZ];
            res2_wavg[a2*NVAR + VAR_DENSITY_ENERGY] = residuals1[a1*NVAR + VAR_DENSITY_ENERGY];
            w_sums[a2] = 1.0;
        } else {
            // Calculate contribution of a1 -> a2:
            const double idist_a1a2 = 1.0/sqrt(dx_a1a2*dx_a1a2 + dy_a1a2*dy_a1a2 + dz_a1a2*dz_a1a2);
            res2_wavg[a2*NVAR + VAR_DENSITY]        += idist_a1a2*residuals1[a1*NVAR + VAR_DENSITY];
            res2_wavg[a2*NVAR + VAR_MOMENTUMX]      += idist_a1a2*residuals1[a1*NVAR + VAR_MOMENTUMX];
            res2_wavg[a2*NVAR + VAR_MOMENTUMY]      += idist_a1a2*residuals1[a1*NVAR + VAR_MOMENTUMY];
            res2_wavg[a2*NVAR + VAR_MOMENTUMZ]      += idist_a1a2*residuals1[a1*NVAR + VAR_MOMENTUMZ];
            res2_wavg[a2*NVAR + VAR_DENSITY_ENERGY] += idist_a1a2*residuals1[a1*NVAR + VAR_DENSITY_ENERGY];
            w_sums[a2] += idist_a1a2;

            // Calculate contribution of b1 >- a2:
            double dx_b1a2 = cb1.x - ca2.x;
            double dy_b1a2 = cb1.y - ca2.y;
            double dz_b1a2 = cb1.z - ca2.z;

            const double idist_b1a2 = 1.0/sqrt(dx_b1a2*dx_b1a2 + dy_b1a2*dy_b1a2 + dz_b1a2*dz_b1a2);
            res2_wavg[a2*NVAR + VAR_DENSITY]        += idist_b1a2*residuals1[b1*NVAR + VAR_DENSITY];
            res2_wavg[a2*NVAR + VAR_MOMENTUMX]      += idist_b1a2*residuals1[b1*NVAR + VAR_MOMENTUMX];
            res2_wavg[a2*NVAR + VAR_MOMENTUMY]      += idist_b1a2*residuals1[b1*NVAR + VAR_MOMENTUMY];
            res2_wavg[a2*NVAR + VAR_MOMENTUMZ]      += idist_b1a2*residuals1[b1*NVAR + VAR_MOMENTUMZ];
            res2_wavg[a2*NVAR + VAR_DENSITY_ENERGY] += idist_b1a2*residuals1[b1*NVAR + VAR_DENSITY_ENERGY];
            w_sums[a2] += idist_b1a2;
        }

        // Process b2:
        double dx_b1b2 = cb2.x - cb1.x;
        double dy_b1b2 = cb2.y - cb1.y;
        double dz_b1b2 = cb2.z - cb1.z;
        if (dx_b1b2 == 0.0 && dy_b1b2 == 0.0 && dz_b1b2 == 0.0) {
            // b2 == b1:
            res2_wavg[b2*NVAR + VAR_DENSITY]        = residuals1[b1*NVAR + VAR_DENSITY];
            res2_wavg[b2*NVAR + VAR_MOMENTUMX]      = residuals1[b1*NVAR + VAR_MOMENTUMX];
            res2_wavg[b2*NVAR + VAR_MOMENTUMY]      = residuals1[b1*NVAR + VAR_MOMENTUMY];
            res2_wavg[b2*NVAR + VAR_MOMENTUMZ]      = residuals1[b1*NVAR + VAR_MOMENTUMZ];
            res2_wavg[b2*NVAR + VAR_DENSITY_ENERGY] = residuals1[b1*NVAR + VAR_DENSITY_ENERGY];
            w_sums[b2] = 1.0;
        } else {
            // Calculate contribution of b1 -> b2:
            const double idist_b1b2 = 1.0/sqrt(dx_b1b2*dx_b1b2 + dy_b1b2*dy_b1b2 + dz_b1b2*dz_b1b2);
            res2_wavg[b2*NVAR + VAR_DENSITY]        += idist_b1b2*residuals1[b1*NVAR + VAR_DENSITY];
            res2_wavg[b2*NVAR + VAR_MOMENTUMX]      += idist_b1b2*residuals1[b1*NVAR + VAR_MOMENTUMX];
            res2_wavg[b2*NVAR + VAR_MOMENTUMY]      += idist_b1b2*residuals1[b1*NVAR + VAR_MOMENTUMY];
            res2_wavg[b2*NVAR + VAR_MOMENTUMZ]      += idist_b1b2*residuals1[b1*NVAR + VAR_MOMENTUMZ];
            res2_wavg[b2*NVAR + VAR_DENSITY_ENERGY] += idist_b1b2*residuals1[b1*NVAR + VAR_DENSITY_ENERGY];
            w_sums[b2] += idist_b1b2;

            // Calculate contribution of a1 -> b2:
            double dx_a1b2 = ca1.x - cb2.x;
            double dy_a1b2 = ca1.y - cb2.y;
            double dz_a1b2 = ca1.z - cb2.z;

            const double idist_a1b2 = 1.0/sqrt(dx_a1b2*dx_a1b2 + dy_a1b2*dy_a1b2 + dz_a1b2*dz_a1b2);
            res2_wavg[b2*NVAR + VAR_DENSITY]        += idist_a1b2*residuals1[b1*NVAR + VAR_DENSITY];
            res2_wavg[b2*NVAR + VAR_MOMENTUMX]      += idist_a1b2*residuals1[b1*NVAR + VAR_MOMENTUMX];
            res2_wavg[b2*NVAR + VAR_MOMENTUMY]      += idist_a1b2*residuals1[b1*NVAR + VAR_MOMENTUMY];
            res2_wavg[b2*NVAR + VAR_MOMENTUMZ]      += idist_a1b2*residuals1[b1*NVAR + VAR_MOMENTUMZ];
            res2_wavg[b2*NVAR + VAR_DENSITY_ENERGY] += idist_a1b2*residuals1[b1*NVAR + VAR_DENSITY_ENERGY];
            w_sums[b2] += idist_a1b2;
        }
    }

    #ifdef TIME
    stop_timer();
    #endif
    #ifdef PAPI
    stop_papi();
    #endif

    #if defined OMP && (defined OMP_SCATTERS)
        }
    #endif



    // 2) Perform the averaging stage, then apply:
    loop_start = 0;
    loop_end = nel2;
    #if defined OMP && (defined OMP_SCATTERS)
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
    record_iters(loop_start, loop_end);

    for (long i=loop_start; i<loop_end; i++) {
        // Divide through by sum of weights:
        for (long j=0; j<NVAR; j++) {
            res2_wavg[i*NVAR +j] /= w_sums[i];

            const long idx = NVAR*i + j;
            variables2[idx] += residuals2[idx] - res2_wavg[idx];
        }
    }

    #ifdef TIME
    stop_timer();
    #endif
    #ifdef PAPI
    stop_papi();
    #endif

    #if defined OMP && (defined OMP_SCATTERS)
        }
    #endif
}
