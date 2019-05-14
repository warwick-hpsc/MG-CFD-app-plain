// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

const int a = edges[i].a;
const int b = edges[i].b;

double ex = edges[i].x;
double ey = edges[i].y;
double ez = edges[i].z;
#ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
    double edge_len = edge_weights[i];
#endif

// Process edge-point A:
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

// Process edge-point B:
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

double p_a_val  = p_b + ex;
double pe_a_val = pe_b + ey;
double mx_a_val = momentum_b.x + ez;
double my_a_val = momentum_b.y;
double mz_a_val = momentum_b.z;
#ifdef FLUX_PRECOMPUTE_EDGE_WEIGHTS
    double p_b_val = p_a + edge_len;
    #else
    double p_b_val = p_a;
#endif
double pe_b_val = pe_a;
double mx_b_val = momentum_a.x;
double my_b_val = momentum_a.y;
double mz_b_val = momentum_a.z;

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

    fluxes[p_a_flx_idx]  += p_a_val;
    fluxes[mx_a_flx_idx] += mx_a_val;
    fluxes[my_a_flx_idx] += my_a_val;
    fluxes[mz_a_flx_idx] += mz_a_val;
    fluxes[pe_a_flx_idx] += pe_a_val;

    fluxes[p_b_flx_idx]  += p_b_val;
    fluxes[mx_b_flx_idx] += mx_b_val;
    fluxes[my_b_flx_idx] += my_b_val;
    fluxes[mz_b_flx_idx] += mz_b_val;
    fluxes[pe_b_flx_idx] += pe_b_val;
#endif
