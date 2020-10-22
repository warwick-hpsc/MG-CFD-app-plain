// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

#ifndef FLUX_WALL_KERNEL_H
#define FLUX_WALL_KERNEL_H

void compute_wall_flux_edge_kernel(
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
    double3 velocity_b, momentum_b, 
        flux_contribution_i_momentum_x_b, 
        flux_contribution_i_momentum_y_b, 
        flux_contribution_i_momentum_z_b,
        flux_contribution_i_density_energy_b;
    p_b          = variables_b[VAR_DENSITY];
    pe_b         = variables_b[VAR_MOMENTUMX];
    momentum_b.x = variables_b[VAR_MOMENTUMY];
    momentum_b.y = variables_b[VAR_MOMENTUMZ];
    momentum_b.z = variables_b[VAR_DENSITY_ENERGY];

    #ifdef FLUX_REUSE_DIV
        compute_velocity_reciprocal(((double)1.0)/p_b, momentum_b, velocity_b);
    #else
        compute_velocity(p_b, momentum_b, velocity_b);
    #endif

    double speed_sqd_b = compute_speed_sqd(velocity_b);
    pressure_b = compute_pressure(p_b, pe_b, speed_sqd_b);
    compute_flux_contribution(momentum_b, pe_b,
        pressure_b, velocity_b, 
        flux_contribution_i_momentum_x_b,
        flux_contribution_i_momentum_y_b,
        flux_contribution_i_momentum_z_b, 
        flux_contribution_i_density_energy_b);

    double factor_x = 0.5*ex;
    double factor_y = 0.5*ey;
    double factor_z = 0.5*ez;

    double p_b_val  = factor_x*(ff_variable[VAR_MOMENTUMX] + momentum_b.x) 
                    + factor_y*(ff_variable[VAR_MOMENTUMY] + momentum_b.y)
                    + factor_z*(ff_variable[VAR_MOMENTUMZ] + momentum_b.z);
    double pe_b_val = factor_x*(ff_flux_contribution_density_energy.x + flux_contribution_i_density_energy_b.x)
                    + factor_y*(ff_flux_contribution_density_energy.y + flux_contribution_i_density_energy_b.y)
                    + factor_z*(ff_flux_contribution_density_energy.z + flux_contribution_i_density_energy_b.z);
    double mx_b_val = factor_x*(ff_flux_contribution_momentum_x.x + flux_contribution_i_momentum_x_b.x)
                    + factor_y*(ff_flux_contribution_momentum_x.y + flux_contribution_i_momentum_x_b.y)
                    + factor_z*(ff_flux_contribution_momentum_x.z + flux_contribution_i_momentum_x_b.z);
    double my_b_val = factor_x*(ff_flux_contribution_momentum_y.x + flux_contribution_i_momentum_y_b.x)
                    + factor_y*(ff_flux_contribution_momentum_y.y + flux_contribution_i_momentum_y_b.y)
                    + factor_z*(ff_flux_contribution_momentum_y.z + flux_contribution_i_momentum_y_b.z);
    double mz_b_val = factor_x*(ff_flux_contribution_momentum_z.x + flux_contribution_i_momentum_z_b.x)
                    + factor_y*(ff_flux_contribution_momentum_z.y + flux_contribution_i_momentum_z_b.y)
                    + factor_z*(ff_flux_contribution_momentum_z.z + flux_contribution_i_momentum_z_b.z);

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