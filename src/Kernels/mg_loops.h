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

#ifndef MG_H
#define MG_H

#include "common.h"

void up(
    double *restrict variables1, 
    double *restrict variables2, 
    int nel2, 
    int *restrict mapping, 
    int *restrict up_scratch, 
    int mgc);

void down(
    double *restrict variables1, 
    int nel1, 
    double *restrict variables2, 
    int nel2, 
    int *restrict mapping, 
    int mgc, 
    double3 *restrict coords1, 
    double3 *restrict coords2);

void down_residuals(
    double *restrict residuals1, 
    int nel1, 
    // double *restrict variables2, 
    // double *restrict residuals2, 
    // Depending on MG configuration, variables2 and residuals2
    // may point to the same array.
    double *variables2, 
    double *residuals2, 
    int nel2, 
    int *restrict mapping, 
    int mgc, 
    double3 *restrict coords1, 
    double3 *restrict coords2);

void down_interpolate(
    double *restrict variables1, 
    int nel1, 
    double *restrict variables2, 
    int nel2, 
    int *restrict mapping, 
    int mgc, 
    double3 *restrict coords1, 
    double3 *restrict coords2);

void down_residuals_interpolate_crude(
    double *restrict residuals1, 
    int nel1, 
    double *restrict residuals2,
    double *restrict variables2, 
    int nel2, 
    int *restrict mapping, 
    int mgc, 
    double3 *restrict coords1, 
    double3 *restrict coords2);

void down_residuals_interpolate_proper(
    edge_neighbour *edges,
    int num_edges,
    double *restrict residuals1, 
    double *restrict residuals2,
    double *restrict variables2, 
    int nel2,
    int *restrict mapping, 
    int mgc, 
    double3 *restrict coords1, 
    double3 *restrict coords2);

#endif
