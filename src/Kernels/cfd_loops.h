// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

#ifndef KERNELS_H
#define KERNELS_H

#ifdef PAPI
#include <papi.h>
#endif

#include "common.h"

void compute_step_factor(
    long nel, 
    const double *restrict variables, 
    const double *restrict volumes, 
    double *restrict step_factors);

void compute_step_factor_legacy(
    long nel, 
    const double *restrict variables, 
    const double *restrict areas, 
    double *restrict step_factors);

void update_edges(
    long first_edge, 
    long nedges, 
    const edge_neighbour *restrict edges,
    const edge *restrict edge_variables,
    double *restrict fluxes);

void time_step(
    int j, 
    long nel, 
    const double *restrict step_factors, 
    double *restrict fluxes, 
    const double *restrict old_variables, 
    double *restrict variables);

void zero_fluxes(
    long nel, 
    double *restrict array);

inline void initialize_variables(
    long nel, 
    double* restrict variables)
{
    log("initialize_variables()");

    for(long i=0; i<nel; i++) {
        for (int v=0; v<NVAR; v++) {
            variables[i*NVAR + v] = ff_variable[v];
        }
    }
}

FORCE_INLINE
inline void compute_flux_contribution(
    double3 momentum, 
    double density_energy, 
    double pressure, 
    double3 velocity, 
    double3& fc_momentum_x, 
    double3& fc_momentum_y, 
    double3& fc_momentum_z, 
    double3& fc_density_energy)
{
    fc_momentum_x.x = velocity.x*momentum.x + pressure;
    fc_momentum_x.y = velocity.x*momentum.y;
    fc_momentum_x.z = velocity.x*momentum.z;

    fc_momentum_y.x = fc_momentum_x.y;
    fc_momentum_y.y = velocity.y*momentum.y + pressure;
    fc_momentum_y.z = velocity.y*momentum.z;

    fc_momentum_z.x = fc_momentum_x.z;
    fc_momentum_z.y = fc_momentum_y.z;
    fc_momentum_z.z = velocity.z*momentum.z + pressure;

    double de_p = density_energy+pressure;
    fc_density_energy.x = velocity.x*de_p;
    fc_density_energy.y = velocity.y*de_p;
    fc_density_energy.z = velocity.z*de_p;
}

inline void initialize_far_field_conditions()
{
    const double angle_of_attack = double(3.1415926535897931 / 180.0) * double(deg_angle_of_attack);

    ff_variable[VAR_DENSITY] = double(1.4);

    double ff_pressure = double(1.0);
    double ff_speed_of_sound = sqrt(GAMMA*ff_pressure / ff_variable[VAR_DENSITY]);
    double ff_speed = double(ff_mach)*ff_speed_of_sound;

    double3 ff_velocity;
    ff_velocity.x = ff_speed*double(cos((double)angle_of_attack));
    ff_velocity.y = ff_speed*double(sin((double)angle_of_attack));
    ff_velocity.z = 0.0;

    ff_variable[VAR_MOMENTUMX] = ff_variable[VAR_DENSITY] * ff_velocity.x;
    ff_variable[VAR_MOMENTUMY] = ff_variable[VAR_DENSITY] * ff_velocity.y;
    ff_variable[VAR_MOMENTUMZ] = ff_variable[VAR_DENSITY] * ff_velocity.z;

    ff_variable[VAR_DENSITY_ENERGY] = ff_variable[VAR_DENSITY]*(double(0.5)*(ff_speed*ff_speed)) + (ff_pressure / double(GAMMA-1.0));

    double3 ff_momentum;
    ff_momentum.x = ff_variable[VAR_MOMENTUMX];
    ff_momentum.y = ff_variable[VAR_MOMENTUMY];
    ff_momentum.z = ff_variable[VAR_MOMENTUMZ];
    compute_flux_contribution(
        ff_momentum, 
        ff_variable[VAR_DENSITY_ENERGY], 
        ff_pressure, 
        ff_velocity, 
        ff_flux_contribution_momentum_x, 
        ff_flux_contribution_momentum_y, 
        ff_flux_contribution_momentum_z, 
        ff_flux_contribution_density_energy);
}

FORCE_INLINE
inline void compute_velocity(double density, double3 momentum, double3& velocity)
{
    velocity.x = momentum.x / density;
    velocity.y = momentum.y / density;
    velocity.z = momentum.z / density;
}

FORCE_INLINE
inline void compute_velocity_reciprocal(double density_reciprocal, double3 momentum, double3& velocity)
{
    velocity.x = momentum.x * density_reciprocal;
    velocity.y = momentum.y * density_reciprocal;
    velocity.z = momentum.z * density_reciprocal;
}

FORCE_INLINE
inline double compute_speed_sqd(double3 velocity)
{
    return velocity.x*velocity.x + velocity.y*velocity.y + velocity.z*velocity.z;
}

FORCE_INLINE
inline double compute_pressure(double density, double density_energy, double speed_sqd)
{
    return (double(GAMMA)-double(1.0))*(density_energy - double(0.5)*density*speed_sqd);
}

FORCE_INLINE
inline double compute_speed_of_sound(double density, double pressure)
{
    return sqrt(double(GAMMA)*pressure/density);
}

FORCE_INLINE
inline double compute_speed_of_sound_reciprocal(double density_reciprocal, double pressure)
{
    return sqrt(double(GAMMA)*pressure*density_reciprocal);
}

FORCE_INLINE
inline void update_a(
    long i,
    long a,
    double *restrict fluxes,
    const edge *restrict edge_variables)
{
    const double p_val  = edge_variables[i*NVAR + VAR_DENSITY].a;
    const double mx_val = edge_variables[i*NVAR + VAR_MOMENTUMX].a;
    const double my_val = edge_variables[i*NVAR + VAR_MOMENTUMY].a;
    const double mz_val = edge_variables[i*NVAR + VAR_MOMENTUMZ].a;
    const double pe_val = edge_variables[i*NVAR + VAR_DENSITY_ENERGY].a;

    const long p_idx  = a*NVAR + VAR_DENSITY;
    const long mx_idx = a*NVAR + VAR_MOMENTUMX;
    const long my_idx = a*NVAR + VAR_MOMENTUMY;
    const long mz_idx = a*NVAR + VAR_MOMENTUMZ;
    const long pe_idx = a*NVAR + VAR_DENSITY_ENERGY;

    fluxes[p_idx]  += p_val;
    fluxes[pe_idx] += pe_val;
    fluxes[mx_idx] += mx_val;
    fluxes[my_idx] += my_val;
    fluxes[mz_idx] += mz_val;
}

FORCE_INLINE
inline void update_b(
    long i,
    long b,
    double *restrict fluxes,
    const edge* restrict edge_variables)
{
    const double p_val  = edge_variables[i*NVAR + VAR_DENSITY].b;
    const double mx_val = edge_variables[i*NVAR + VAR_MOMENTUMX].b;
    const double my_val = edge_variables[i*NVAR + VAR_MOMENTUMY].b;
    const double mz_val = edge_variables[i*NVAR + VAR_MOMENTUMZ].b;
    const double pe_val = edge_variables[i*NVAR + VAR_DENSITY_ENERGY].b;

    const long p_idx  = b*NVAR + VAR_DENSITY;
    const long mx_idx = b*NVAR + VAR_MOMENTUMX;
    const long my_idx = b*NVAR + VAR_MOMENTUMY;
    const long mz_idx = b*NVAR + VAR_MOMENTUMZ;
    const long pe_idx = b*NVAR + VAR_DENSITY_ENERGY;

    fluxes[p_idx]  += p_val;
    fluxes[pe_idx] += pe_val;
    fluxes[mx_idx] += mx_val;
    fluxes[my_idx] += my_val;
    fluxes[mz_idx] += mz_val;
}

#endif
