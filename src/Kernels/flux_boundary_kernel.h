// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

#ifndef COMPUTE_BOUNDARY_FLUX_EDGE_KERNEL_H
#define COMPUTE_BOUNDARY_FLUX_EDGE_KERNEL_H

void compute_boundary_flux_edge_kernel(
    double ex, double ey, double ez,
    const double *restrict variables_b, 
    #ifdef FLUX_FISSION
        edge *restrict edge_variables
    #else
        double *restrict fluxes_b
    #endif
    )
{
    double p_b, pe_b, pressure_b;
    double3 velocity_b, momentum_b;
    p_b          = variables_b[VAR_DENSITY];
    pe_b         = variables_b[VAR_MOMENTUMX];
    momentum_b.x = variables_b[VAR_MOMENTUMY];
    momentum_b.y = variables_b[VAR_MOMENTUMZ];
    momentum_b.z = variables_b[VAR_DENSITY_ENERGY];

    compute_velocity(p_b, momentum_b, velocity_b);

    double speed_sqd_b = compute_speed_sqd(velocity_b);
    pressure_b = compute_pressure(p_b, pe_b, speed_sqd_b);

    double p_b_val = 0;
    double mx_b_val = ex*pressure_b;
    double my_b_val = ey*pressure_b;
    double mz_b_val = ez*pressure_b;
    double pe_b_val = 0;

    #ifdef FLUX_FISSION
        edge_variables[VAR_DENSITY].b   = p_b_val;
        edge_variables[VAR_MOMENTUMX].b = mx_b_val;
        edge_variables[VAR_MOMENTUMY].b = my_b_val;
        edge_variables[VAR_MOMENTUMZ].b = mz_b_val;
        edge_variables[VAR_DENSITY_ENERGY].b = pe_b_val;
    #else
        fluxes_b[VAR_DENSITY]  += p_b_val;
        fluxes_b[VAR_MOMENTUMX] += mx_b_val;
        fluxes_b[VAR_MOMENTUMY] += my_b_val;
        fluxes_b[VAR_MOMENTUMZ] += mz_b_val;
        fluxes_b[VAR_DENSITY_ENERGY] += pe_b_val;
    #endif
}

#endif
