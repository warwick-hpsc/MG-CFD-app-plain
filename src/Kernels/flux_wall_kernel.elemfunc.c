// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

// void compute_wall_flux_edge(
//     long first_edge,
//     long nedges,
//     const edge_neighbour *restrict edges, 
//     const double *restrict variables, 
//     #ifdef FLUX_FISSION
//         edge *restrict edge_variables
//     #else
//         double *restrict fluxes
//     #endif)

long b = edges[i].b;

double p_b, pe_b, pressure_b;
double3 velocity_b, momentum_b, 
    flux_contribution_i_momentum_x_b, 
    flux_contribution_i_momentum_y_b, 
    flux_contribution_i_momentum_z_b,
    flux_contribution_i_density_energy_b;

const long p_b_idx  = b*NVAR + VAR_DENSITY;
const long mx_b_idx = b*NVAR + VAR_MOMENTUMX;
const long my_b_idx = b*NVAR + VAR_MOMENTUMY;
const long mz_b_idx = b*NVAR + VAR_MOMENTUMZ;
const long pe_b_idx = b*NVAR + VAR_DENSITY_ENERGY;

p_b  = variables[p_b_idx];
pe_b = variables[pe_b_idx];
momentum_b.x = variables[mx_b_idx];
momentum_b.y = variables[my_b_idx];
momentum_b.z = variables[mz_b_idx];

#ifdef FLUX_REUSE_DIV
    compute_velocity_reciprocal(((double)1.0)/p_b, momentum_b, velocity_b);
#else
    compute_velocity(p_b, momentum_b, velocity_b);
#endif

double speed_sqd_b = compute_speed_sqd(velocity_b);
pressure_b = compute_pressure(p_b, pe_b, speed_sqd_b);
compute_flux_contribution(p_b, momentum_b, pe_b,
    pressure_b, velocity_b, 
    flux_contribution_i_momentum_x_b,
    flux_contribution_i_momentum_y_b,
    flux_contribution_i_momentum_z_b, 
    flux_contribution_i_density_energy_b);

double factor_x = 0.5*edges[i].x;
double factor_y = 0.5*edges[i].y;
double factor_z = 0.5*edges[i].z;

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
    edge_variables[i*NVAR + VAR_DENSITY].b   = p_b_val;
    edge_variables[i*NVAR + VAR_MOMENTUMX].b = mx_b_val;
    edge_variables[i*NVAR + VAR_MOMENTUMY].b = my_b_val;
    edge_variables[i*NVAR + VAR_MOMENTUMZ].b = mz_b_val;
    edge_variables[i*NVAR + VAR_DENSITY_ENERGY].b = pe_b_val;
#else
    const long p_b_flx_idx  = b*NVAR + VAR_DENSITY;
    const long mx_b_flx_idx = b*NVAR + VAR_MOMENTUMX;
    const long my_b_flx_idx = b*NVAR + VAR_MOMENTUMY;
    const long mz_b_flx_idx = b*NVAR + VAR_MOMENTUMZ;
    const long pe_b_flx_idx = b*NVAR + VAR_DENSITY_ENERGY;
    
    fluxes[p_b_flx_idx]  += p_b_val;
    fluxes[mx_b_flx_idx] += mx_b_val;
    fluxes[my_b_flx_idx] += my_b_val;
    fluxes[mz_b_flx_idx] += mz_b_val;
    fluxes[pe_b_flx_idx] += pe_b_val;
#endif
