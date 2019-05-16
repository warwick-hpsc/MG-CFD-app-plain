// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

// void compute_flux_edge(
//     int first_edge,
//     int nedges,
//     const edge_neighbour *restrict edges, 
//     #ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
//         const double *restrict edge_weights,
//     #endif
//     const double *restrict variables, 
//     #ifdef FLUX_FISSION
//         edge *restrict edge_variables
//     #else
//         double *restrict fluxes
//     #endif)

const int a = edges[i].a;
const int b = edges[i].b;

double ex = edges[i].x;
double ey = edges[i].y;
double ez = edges[i].z;
#ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
    double ewt = edge_weights[i];
#else
    double ewt = sqrt(ex*ex + ey*ey + ez*ez);
#endif

////////////////////////////////////
// Read and process edge-point B:
////////////////////////////////////
double p_b, pe_b;
double3 momentum_b;
const int p_b_idx  = b*NVAR + VAR_DENSITY;
const int mx_b_idx = b*NVAR + VAR_MOMENTUMX;
const int my_b_idx = b*NVAR + VAR_MOMENTUMY;
const int mz_b_idx = b*NVAR + VAR_MOMENTUMZ;
const int pe_b_idx = b*NVAR + VAR_DENSITY_ENERGY;
p_b          = variables[ p_b_idx];
momentum_b.x = variables[mx_b_idx];
momentum_b.y = variables[my_b_idx];
momentum_b.z = variables[mz_b_idx];
pe_b         = variables[pe_b_idx];

#ifdef FLUX_REUSE_DIV
double p_b_reciprocal = ((double)1.0)/p_b;
#endif

double3 velocity_b,
    flux_contribution_i_momentum_x_b, 
    flux_contribution_i_momentum_y_b, 
    flux_contribution_i_momentum_z_b,
    flux_contribution_i_density_energy_b;

#ifdef FLUX_REUSE_DIV
    // Save 3x divisions:
    compute_velocity_reciprocal(p_b_reciprocal, momentum_b, velocity_b);
#else
    compute_velocity(p_b, momentum_b, velocity_b);
#endif
double speed_sqd_b = compute_speed_sqd(velocity_b);
double speed_b = sqrt(speed_sqd_b);

double pressure_b = compute_pressure(p_b, pe_b, speed_sqd_b);
#ifdef FLUX_REUSE_DIV
    // Save 1x division:
    double speed_of_sound_b = compute_speed_of_sound_reciprocal(p_b_reciprocal, pressure_b);
#else
    double speed_of_sound_b = compute_speed_of_sound(p_b, pressure_b);
#endif
compute_flux_contribution(p_b, momentum_b, pe_b,
                          pressure_b, velocity_b, 
                          flux_contribution_i_momentum_x_b,
                          flux_contribution_i_momentum_y_b,
                          flux_contribution_i_momentum_z_b, 
                          flux_contribution_i_density_energy_b);

////////////////////////////////////
// Read and process edge-point A:
////////////////////////////////////
double p_a, pe_a;
double3 momentum_a;
const int p_a_idx  = a*NVAR + VAR_DENSITY;
const int mx_a_idx = a*NVAR + VAR_MOMENTUMX;
const int my_a_idx = a*NVAR + VAR_MOMENTUMY;
const int mz_a_idx = a*NVAR + VAR_MOMENTUMZ;
const int pe_a_idx = a*NVAR + VAR_DENSITY_ENERGY;
p_a          = variables[ p_a_idx];
momentum_a.x = variables[mx_a_idx];
momentum_a.y = variables[my_a_idx];
momentum_a.z = variables[mz_a_idx];
pe_a         = variables[pe_a_idx];

#ifdef FLUX_REUSE_DIV
    double p_a_reciprocal = ((double)1.0)/p_a;
#endif

double3 velocity_a, 
    flux_contribution_i_momentum_x_a, 
    flux_contribution_i_momentum_y_a, 
    flux_contribution_i_momentum_z_a,
    flux_contribution_i_density_energy_a;
#ifdef FLUX_REUSE_DIV
    // Save 3x divisions:
    compute_velocity_reciprocal(p_a_reciprocal, momentum_a, velocity_a);
#else
    compute_velocity(p_a, momentum_a, velocity_a);
#endif

double speed_sqd_a = compute_speed_sqd(velocity_a);
double speed_a = sqrt(speed_sqd_a);

double pressure_a = compute_pressure(p_a, pe_a, speed_sqd_a);

#ifdef FLUX_REUSE_DIV
    // Save 1x division:
    double speed_of_sound_a = compute_speed_of_sound_reciprocal(p_a_reciprocal, pressure_a);
#else
    double speed_of_sound_a = compute_speed_of_sound(p_a, pressure_a);
#endif

compute_flux_contribution(p_a, momentum_a, pe_a,
                        pressure_a, velocity_a, 
                        flux_contribution_i_momentum_x_a,
                        flux_contribution_i_momentum_y_a,
                        flux_contribution_i_momentum_z_a, 
                        flux_contribution_i_density_energy_a);

double factor_a = -ewt*smoothing_coefficient*0.5
                  *(speed_a + speed_b + speed_of_sound_a + speed_of_sound_b);
#ifdef FLUX_REUSE_FACTOR
    double factor_b = factor_a;
#else
    double factor_b = -ewt*smoothing_coefficient*0.5
                    *(speed_a + speed_b + speed_of_sound_a + speed_of_sound_b);
#endif

double factor_x = -0.5*ex;
double factor_y = -0.5*ey;
double factor_z = -0.5*ez;

double p_a_val  = factor_a*(p_a - p_b)
                + factor_x*(momentum_a.x + momentum_b.x)       
                + factor_y*(momentum_a.y + momentum_b.y)
                + factor_z*(momentum_a.z + momentum_b.z);
double pe_a_val = factor_a*(pe_a - pe_b)
                + factor_x*(flux_contribution_i_density_energy_a.x + flux_contribution_i_density_energy_b.x)
                + factor_y*(flux_contribution_i_density_energy_a.y + flux_contribution_i_density_energy_b.y)
                + factor_z*(flux_contribution_i_density_energy_a.z + flux_contribution_i_density_energy_b.z);
double mx_a_val = factor_a*(momentum_a.x - momentum_b.x) 
                + factor_x*(flux_contribution_i_momentum_x_a.x + flux_contribution_i_momentum_x_b.x) 
                + factor_y*(flux_contribution_i_momentum_x_a.y + flux_contribution_i_momentum_x_b.y)
                + factor_z*(flux_contribution_i_momentum_x_a.z + flux_contribution_i_momentum_x_b.z);
double my_a_val = factor_a*(momentum_a.y - momentum_b.y)
                + factor_x*(flux_contribution_i_momentum_y_a.x + flux_contribution_i_momentum_y_b.x)
                + factor_y*(flux_contribution_i_momentum_y_a.y + flux_contribution_i_momentum_y_b.y)
                + factor_z*(flux_contribution_i_momentum_y_a.z + flux_contribution_i_momentum_y_b.z);
double mz_a_val = factor_a*(momentum_a.z - momentum_b.z)
                + factor_x*(flux_contribution_i_momentum_z_a.x + flux_contribution_i_momentum_z_b.x)
                + factor_y*(flux_contribution_i_momentum_z_a.y + flux_contribution_i_momentum_z_b.y)
                + factor_z*(flux_contribution_i_momentum_z_a.z + flux_contribution_i_momentum_z_b.z);
#ifdef FLUX_REUSE_FLUX
    double p_b_val  = -p_a_val;
    double pe_b_val = -pe_a_val;
    double mx_b_val = -mx_a_val;
    double my_b_val = -my_a_val;
    double mz_b_val = -mz_a_val;
#else
    double p_b_val  = factor_b*(p_b - p_a)
                    - factor_x*(momentum_a.x + momentum_b.x)
                    - factor_y*(momentum_a.y + momentum_b.y)
                    - factor_z*(momentum_a.z + momentum_b.z);
    double pe_b_val = factor_b*(pe_b - pe_a)
                    - factor_x*(flux_contribution_i_density_energy_a.x + flux_contribution_i_density_energy_b.x)
                    - factor_y*(flux_contribution_i_density_energy_a.y + flux_contribution_i_density_energy_b.y)
                    - factor_z*(flux_contribution_i_density_energy_a.z + flux_contribution_i_density_energy_b.z);
    double mx_b_val = factor_b*(momentum_b.x - momentum_a.x)
                    - factor_x*(flux_contribution_i_momentum_x_a.x + flux_contribution_i_momentum_x_b.x)
                    - factor_y*(flux_contribution_i_momentum_x_a.y + flux_contribution_i_momentum_x_b.y)
                    - factor_z*(flux_contribution_i_momentum_x_a.z + flux_contribution_i_momentum_x_b.z);
    double my_b_val = factor_b*(momentum_b.y - momentum_a.y)    
                    - factor_x*(flux_contribution_i_momentum_y_a.x + flux_contribution_i_momentum_y_b.x)
                    - factor_y*(flux_contribution_i_momentum_y_a.y + flux_contribution_i_momentum_y_b.y)
                    - factor_z*(flux_contribution_i_momentum_y_a.z + flux_contribution_i_momentum_y_b.z);
    double mz_b_val = factor_b*(momentum_b.z - momentum_a.z)
                    - factor_x*(flux_contribution_i_momentum_z_a.x + flux_contribution_i_momentum_z_b.x)
                    - factor_y*(flux_contribution_i_momentum_z_a.y + flux_contribution_i_momentum_z_b.y)
                    - factor_z*(flux_contribution_i_momentum_z_a.z + flux_contribution_i_momentum_z_b.z);
#endif

// Write out fluxes:
#ifdef FLUX_FISSION
    edge_variables[i*NVAR + VAR_DENSITY       ].a =  p_a_val;
    edge_variables[i*NVAR + VAR_MOMENTUMX     ].a = mx_a_val;
    edge_variables[i*NVAR + VAR_MOMENTUMY     ].a = my_a_val;
    edge_variables[i*NVAR + VAR_MOMENTUMZ     ].a = mz_a_val;
    edge_variables[i*NVAR + VAR_DENSITY_ENERGY].a = pe_a_val;

    edge_variables[i*NVAR + VAR_DENSITY       ].b =  p_b_val;
    edge_variables[i*NVAR + VAR_MOMENTUMX     ].b = mx_b_val;
    edge_variables[i*NVAR + VAR_MOMENTUMY     ].b = my_b_val;
    edge_variables[i*NVAR + VAR_MOMENTUMZ     ].b = mz_b_val;
    edge_variables[i*NVAR + VAR_DENSITY_ENERGY].b = pe_b_val;
#else
    const int p_a_flx_idx  = a*NVAR + VAR_DENSITY;
    const int mx_a_flx_idx = a*NVAR + VAR_MOMENTUMX;
    const int my_a_flx_idx = a*NVAR + VAR_MOMENTUMY;
    const int mz_a_flx_idx = a*NVAR + VAR_MOMENTUMZ;
    const int pe_a_flx_idx = a*NVAR + VAR_DENSITY_ENERGY;

    const int p_b_flx_idx  = b*NVAR + VAR_DENSITY;
    const int mx_b_flx_idx = b*NVAR + VAR_MOMENTUMX;
    const int my_b_flx_idx = b*NVAR + VAR_MOMENTUMY;
    const int mz_b_flx_idx = b*NVAR + VAR_MOMENTUMZ;
    const int pe_b_flx_idx = b*NVAR + VAR_DENSITY_ENERGY;

    fluxes[p_a_flx_idx]  +=  p_a_val;
    fluxes[mx_a_flx_idx] += mx_a_val;
    fluxes[my_a_flx_idx] += my_a_val;
    fluxes[mz_a_flx_idx] += mz_a_val;
    fluxes[pe_a_flx_idx] += pe_a_val;

    fluxes[p_b_flx_idx]  +=  p_b_val;
    fluxes[mx_b_flx_idx] += mx_b_val;
    fluxes[my_b_flx_idx] += my_b_val;
    fluxes[mz_b_flx_idx] += mz_b_val;
    fluxes[pe_b_flx_idx] += pe_b_val;
#endif
