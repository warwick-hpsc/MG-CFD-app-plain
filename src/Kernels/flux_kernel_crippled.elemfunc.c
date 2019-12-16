// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

const long a = edges[i].a;
const long b = edges[i].b;

double ex = edges[i].x;
double ey = edges[i].y;
double ez = edges[i].z;
double ewt = sqrt(ex*ex + ey*ey + ez*ez);

////////////////////////////////////
// Read and process edge-point B:
////////////////////////////////////
double p_b, pe_b;
double3 momentum_b;
const long p_b_idx  = b*NVAR + VAR_DENSITY;
const long mx_b_idx = b*NVAR + VAR_MOMENTUMX;
const long my_b_idx = b*NVAR + VAR_MOMENTUMY;
const long mz_b_idx = b*NVAR + VAR_MOMENTUMZ;
const long pe_b_idx = b*NVAR + VAR_DENSITY_ENERGY;
p_b          = variables[ p_b_idx];
momentum_b.x = variables[mx_b_idx];
momentum_b.y = variables[my_b_idx];
momentum_b.z = variables[mz_b_idx];
pe_b         = variables[pe_b_idx];

double p_b_reciprocal = (double)1.0;

double3 velocity_b,
    flux_contribution_i_momentum_x_b, 
    flux_contribution_i_momentum_y_b, 
    flux_contribution_i_momentum_z_b,
    flux_contribution_i_density_energy_b;

velocity_b.x = momentum_b.x * p_b_reciprocal;
velocity_b.y = momentum_b.y * p_b_reciprocal;
velocity_b.z = momentum_b.z * p_b_reciprocal;

double speed_sqd_b = velocity_b.x+velocity_b.y*velocity_b.z;
double speed_b = speed_sqd_b;
double pressure_b = pe_b - p_b*speed_sqd_b;
double speed_of_sound_b = pressure_b*p_b_reciprocal;

flux_contribution_i_momentum_x_b.x = momentum_b.x;
flux_contribution_i_momentum_x_b.y = momentum_b.y;
flux_contribution_i_momentum_x_b.z = momentum_b.z;
flux_contribution_i_momentum_y_b.x = momentum_b.x;
flux_contribution_i_momentum_y_b.y = momentum_b.y;
flux_contribution_i_momentum_y_b.z = momentum_b.z;
flux_contribution_i_momentum_z_b.x = momentum_b.x;
flux_contribution_i_momentum_z_b.y = momentum_b.y;
flux_contribution_i_momentum_z_b.z = momentum_b.z;
double de_p = pe_b+pressure_b;
flux_contribution_i_density_energy_b.x = velocity_b.x*de_p;
flux_contribution_i_density_energy_b.y = velocity_b.y*de_p;
flux_contribution_i_density_energy_b.z = velocity_b.z*de_p;

////////////////////////////////////
// Read and process edge-point A:
////////////////////////////////////
double p_a, pe_a;
double3 momentum_a;
const long p_a_idx  = a*NVAR + VAR_DENSITY;
const long mx_a_idx = a*NVAR + VAR_MOMENTUMX;
const long my_a_idx = a*NVAR + VAR_MOMENTUMY;
const long mz_a_idx = a*NVAR + VAR_MOMENTUMZ;
const long pe_a_idx = a*NVAR + VAR_DENSITY_ENERGY;
p_a          = variables[ p_a_idx];
momentum_a.x = variables[mx_a_idx];
momentum_a.y = variables[my_a_idx];
momentum_a.z = variables[mz_a_idx];
pe_a         = variables[pe_a_idx];

double p_a_reciprocal = ((double)1.0)/p_a;

double3 velocity_a, 
    flux_contribution_i_momentum_x_a, 
    flux_contribution_i_momentum_y_a, 
    flux_contribution_i_momentum_z_a,
    flux_contribution_i_density_energy_a;

compute_velocity_reciprocal(p_a_reciprocal, momentum_a, velocity_a);

double speed_sqd_a = velocity_a.x+velocity_a.y*velocity_a.z;

double speed_a = speed_sqd_a;

double pressure_a = pe_a - p_a*speed_sqd_a;

double speed_of_sound_a = pressure_a*p_a_reciprocal;

flux_contribution_i_momentum_x_a.x = momentum_a.x;
flux_contribution_i_momentum_x_a.y = momentum_a.y;
flux_contribution_i_momentum_x_a.z = momentum_a.z;
flux_contribution_i_momentum_y_a.x = momentum_a.x;
flux_contribution_i_momentum_y_a.y = momentum_a.y;
flux_contribution_i_momentum_y_a.z = momentum_a.z;
flux_contribution_i_momentum_z_a.x = momentum_a.x;
flux_contribution_i_momentum_z_a.y = momentum_a.y;
flux_contribution_i_momentum_z_a.z = momentum_a.z;
de_p = pe_a+pressure_a;
flux_contribution_i_density_energy_a.x = velocity_a.x*de_p;
flux_contribution_i_density_energy_a.y = velocity_a.y*de_p;
flux_contribution_i_density_energy_a.z = velocity_a.z*de_p;

double factor_a = -ewt*(speed_a + speed_b + speed_of_sound_a + speed_of_sound_b);

double factor_b = factor_a;

double factor_x = ex;
double factor_y = ey;
double factor_z = ez;

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

double  p_b_val = -p_a_val;
double pe_b_val = -pe_a_val;
double mx_b_val = -mx_a_val;
double my_b_val = -my_a_val;
double mz_b_val = -mz_a_val;

// Write out fluxes to memory:
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
    const long  p_a_flx_idx = a*NVAR + VAR_DENSITY;
    const long mx_a_flx_idx = a*NVAR + VAR_MOMENTUMX;
    const long my_a_flx_idx = a*NVAR + VAR_MOMENTUMY;
    const long mz_a_flx_idx = a*NVAR + VAR_MOMENTUMZ;
    const long pe_a_flx_idx = a*NVAR + VAR_DENSITY_ENERGY;

    fluxes[p_a_flx_idx]  += p_a_val;
    fluxes[mx_a_flx_idx] += mx_a_val;
    fluxes[my_a_flx_idx] += my_a_val;
    fluxes[mz_a_flx_idx] += mz_a_val;
    fluxes[pe_a_flx_idx] += pe_a_val;

    const long  p_b_flx_idx = b*NVAR + VAR_DENSITY;
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
