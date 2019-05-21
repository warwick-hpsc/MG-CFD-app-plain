// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

// void compute_boundary_flux_edge(
//     int first_edge,
//     int nedges, 
//     const edge_neighbour *restrict edges, 
//     const double *restrict variables, 
//     #ifdef FLUX_FISSION
//         edge *restrict edge_variables
//     #else
//         double *restrict fluxes
//     #endif)

int b = edges[i].b;

double p_b, pe_b, pressure_b;
double3 velocity_b, momentum_b;

const int p_b_idx  = b*NVAR + VAR_DENSITY;
const int mx_b_idx = b*NVAR + VAR_MOMENTUMX;
const int my_b_idx = b*NVAR + VAR_MOMENTUMY;
const int mz_b_idx = b*NVAR + VAR_MOMENTUMZ;
const int pe_b_idx = b*NVAR + VAR_DENSITY_ENERGY;

p_b  = variables[p_b_idx];
pe_b = variables[pe_b_idx];
momentum_b.x = variables[mx_b_idx];
momentum_b.y = variables[my_b_idx];
momentum_b.z = variables[mz_b_idx];

compute_velocity(p_b, momentum_b, velocity_b);

double speed_sqd_b = compute_speed_sqd(velocity_b);
pressure_b = compute_pressure(p_b, pe_b, speed_sqd_b);

double x = edges[i].x;
double y = edges[i].y;
double z = edges[i].z;

double p_b_val = 0;
double mx_b_val = x*pressure_b;
double my_b_val = y*pressure_b;
double mz_b_val = z*pressure_b;
double pe_b_val = 0;

#ifdef FLUX_FISSION
    edge_variables[i*NVAR + VAR_DENSITY].b   = p_b_val;
    edge_variables[i*NVAR + VAR_MOMENTUMX].b = mx_b_val;
    edge_variables[i*NVAR + VAR_MOMENTUMY].b = my_b_val;
    edge_variables[i*NVAR + VAR_MOMENTUMZ].b = mz_b_val;
    edge_variables[i*NVAR + VAR_DENSITY_ENERGY].b = pe_b_val;
#else
    const int p_b_flx_idx  = b*NVAR + VAR_DENSITY;
    const int mx_b_flx_idx = b*NVAR + VAR_MOMENTUMX;
    const int my_b_flx_idx = b*NVAR + VAR_MOMENTUMY;
    const int mz_b_flx_idx = b*NVAR + VAR_MOMENTUMZ;
    const int pe_b_flx_idx = b*NVAR + VAR_DENSITY_ENERGY;
    
    fluxes[p_b_flx_idx]  += p_b_val;
    fluxes[mx_b_flx_idx] += mx_b_val;
    fluxes[my_b_flx_idx] += my_b_val;
    fluxes[mz_b_flx_idx] += mz_b_val;
    fluxes[pe_b_flx_idx] += pe_b_val;
#endif
